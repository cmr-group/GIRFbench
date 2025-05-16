import numpy as np
import numpy.fft as fft

# TODO: my naming of 'x' and 'y' is completely inconsistent between different functions

class NumericalGIRF:
    def __init__(self, all_x, all_y, extra_op = (), extra_y = (), weight=None, extra_scale = 10, lam = 0):
        self.extra_scale = extra_scale
        
        self.all_x = all_x
        self.all_y = all_y
        
        self.lam = lam
        
        if self.all_x.shape[1] != all_y.shape[1]:
            print('all_x.shape[1] != all_y.shape[1]')
        
        self.all_conv = []
        for x in all_x:
            self.all_conv.append(Convolver(x))
        for op in extra_op:
            self.all_conv.append(op)
         
        self.Ax_size = sum([_.A_shape[0] for _ in self.all_conv])
        self.Ax = np.zeros(self.Ax_size, dtype=all_y.dtype)
        
        # TODO: confirm this matches all op.A_shape[1]
        self.Atx = np.zeros(all_x.shape[1])
        
        if weight is None:
            self.weight = np.ones(all_x.shape[1])
        else:
            self.weight = weight
            
        self.all_y = np.hstack([all_y.ravel(), *extra_y])
        
        if self.Ax_size != self.all_y.size:
            print('WARNING: Ax_size != all_y size.', self.Ax_size, all_y.size)
            
        self.get_spec_norms_normal()
        for i in range(len(extra_op)):
            self.spec_norms[all_x.shape[0]+i] /= extra_scale
            
        
            
    def matvec(self, x):
        # Computes Ax
        start_idx = 0
        stop_idx = 0
        for conv in self.all_conv:
            stop_idx += conv.A_shape[0]
            self.Ax[start_idx:stop_idx] = conv.matvec(x)
            start_idx = stop_idx
            
        return self.Ax
    
    def rmatvec(self, x, spec_norm = False):
        # Computes Ax
        self.Atx *= 0
        
        start_idx = 0
        stop_idx = 0
        for idx, conv in enumerate(self.all_conv):
            stop_idx += conv.A_shape[0]
            if spec_norm:
                self.Atx += conv.rmatvec(x[start_idx:stop_idx]) / self.spec_norms[idx]
            else:
                self.Atx += conv.rmatvec(x[start_idx:stop_idx])
            start_idx = stop_idx
            
        return self.Atx
    
    def AtA(self, x, do_weight = True):
        if do_weight:
            return self.weight*(self.rmatvec(self.matvec(x), spec_norm=True) + self.lam)
        return self.rmatvec(self.matvec(x), spec_norm=True) + self.lam
    
    def Atb(self, x, do_weight = True):
        if do_weight:
            return self.weight*self.rmatvec(x, spec_norm=True)
        return self.rmatvec(x, spec_norm=True)
    
    def get_spec_norms_normal(self):
        
        np.random.seed(42)
        
        self.spec_norms = []
        for conv in self.all_conv:
            x = np.random.randn(self.Atx.size)
            for _ in range(20):
                x_1 = self.weight*conv.rmatvec(conv.matvec(x))
                x_1_norm = np.linalg.norm(x_1)
                x = x_1 / x_1_norm
            self.spec_norms.append(x_1_norm)
            
            
class BGPCResp:
    def __init__(self, g_out, Nsamp, TE_idx, excite_idx, shift):
        self.Nsamp = Nsamp
        self.g_out = g_out
        
        self.N_TR = g_out.size//2

        self.Npad = int(self.Nsamp)
        self.NsampPad = self.Nsamp + 2*self.Npad
        
        self.gg_tile = np.tile(g_out, int(np.ceil(self.NsampPad/g_out.size)))
        self.gg_tile = self.gg_tile[-self.NsampPad:]
        
        self.TE_idx = TE_idx
        self.win_width = TE_idx.max()-TE_idx.min()

        self.shift = shift
        
        self.excite_idx = excite_idx
        # TODO: Calculate all complete 2*TR segments in the data, given the shifting
        # involved
        # We always line up the end og the gradient waveforms with the end of the output
        self.i_excite0 = self.Nsamp-4*self.N_TR + excite_idx + self.shift
        self.i_excite1 = self.i_excite0 + self.N_TR
        
        
        self.win_min = self.TE_idx.min() - self.excite_idx
        self.win_max = self.TE_idx.max() - self.excite_idx
        
        self.possible_windows = []
        for i in range(10):
            p_window = [-i*2*self.N_TR + excite_idx + self.shift, 0]
            p_window[1] = int(p_window[0] + self.N_TR + self.win_max)
            if p_window[1] < 0 and p_window[0] > -self.Nsamp:
                self.possible_windows.append(p_window)
        self.possible_windows = np.array(self.possible_windows)
        
        # print(f'{self.i_excite0 = }')
        # print(f'{self.i_excite1 = }')
        # print(f'{self.possible_windows = }')
        self.possible_windows += self.Nsamp
        # print(f'{self.possible_windows = }')
        
        self.i_excite0 = self.possible_windows[0,0]
        self.i_excite1 = self.i_excite0 + self.N_TR
        
        # print(f'{self.i_excite0 = }')
        # print(f'{self.i_excite1 = }')
        
        
        self.A_shape = (self.win_width, self.Nsamp)

        self.INPUT = fft.fft(self.gg_tile)
        self.INPUTc = np.conj(self.INPUT)
        self.xpad = np.zeros(self.NsampPad, self.gg_tile.dtype)
        self._x = np.zeros(self.Nsamp, self.gg_tile.dtype)
        
        
    def matvec(self, x):
        self.xpad *= 0
        self.xpad[:-2*self.Npad] = x
        y = np.real(fft.ifft(self.INPUT*fft.fft(self.xpad)))[2*self.Npad:]
        
        c_ph0 = np.cumsum(y[self.i_excite0:self.i_excite0+self.win_max])
        c_ph1 = np.cumsum(y[self.i_excite1:self.i_excite1+self.win_max])
        
        resp = c_ph1 - c_ph0
        
        return resp[self.win_min:]
    
    def convolve_only(self, x):
        self.xpad[:-2*self.Npad] = x
        y = np.real(fft.ifft(self.INPUT*fft.fft(self.xpad)))[2*self.Npad:]
        
        return y
    
    def s_transpose(self, resp):
        
        self._x *= 0
        
        c_ph = 0
        for i in range(self.win_max):
            if i < self.win_width:
                c_ph += resp[self.win_width-i-1]

            self._x[self.i_excite0+self.win_max-i] = -c_ph
            self._x[self.i_excite1+self.win_max-i] = c_ph
        
        return self._x

    def rmatvec(self, resp):
        
        self._x *= 0
        
        c_ph = 0
        for i in range(self.win_max):
            if i < self.win_width:
                c_ph += resp[self.win_width-i-1]

            self._x[self.i_excite0+self.win_max-i] = -c_ph
            self._x[self.i_excite1+self.win_max-i] = c_ph
        
        self.xpad *= 0
        self.xpad[2*self.Npad:] = self._x
        return np.real(fft.ifft(self.INPUTc*fft.fft(self.xpad)))[:-2*self.Npad]
        

class Convolver:
    def __init__(self, input, mode='fftconvolve'):
        self.Nsamp = input.size
        self.input = input
        self.mode = mode
        self.A_shape = (self.Nsamp, self.Nsamp)
        
        if self.mode == 'matrix':
            self.A = np.zeros([self.Nsamp, self.Nsamp], input.dtype)
            for i in range(1, self.Nsamp+1):
                self.A[i-1, :i] = input[i-1::-1]
        elif self.mode == 'npconvolve':
            self.inputc = np.real(fft.ifft(np.conj(fft.fft(input))))
            self.Atx = np.zeros_like(input)
        elif self.mode == 'fftconvolve':
            self.Npad = int(self.Nsamp//2)
            self.INPUT = fft.fft(np.pad(input,self.Npad))
            self.INPUTc = np.conj(self.INPUT)
            self.xpad = np.zeros(input.size + 2*self.Npad, input.dtype)
        
        
    def matvec(self, x):
        if self.mode == 'matrix':
            return self.A@x
        if self.mode == 'npconvolve':
            return np.convolve(self.input, x)[:self.Nsamp]
        if self.mode == 'fftconvolve':
            self.xpad[self.Npad:-self.Npad] = x
            return np.real(fft.ifft(self.INPUT*fft.fft(self.xpad)))[2*self.Npad:]
        
        return 0
    
    def rmatvec(self, x):
        if self.mode == 'matrix':
            return self.A.T@x
        if self.mode == 'npconvolve':
            self.Atx[:-1] = np.convolve(self.inputc, x)[self.Nsamp:]
            return self.Atx
        if self.mode == 'fftconvolve':
            self.xpad[self.Npad:-self.Npad] = x
            return np.real(fft.ifft(self.INPUTc*fft.fft(self.xpad)))[:2*self.Npad]
        
        return 0
        
        
        
        
        
        
        
# if self.mode == 'matrix':
#     pass
# elif self.mode == 'npconvolve':
#     pass
# elif self.mode == 'fftconvolve':
#     pass
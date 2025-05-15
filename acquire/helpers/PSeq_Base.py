import datetime

import numpy as np
import pypulseq as pp
from pypulseq.opts import Opts

class PSeq_Params:
    """Hold paramaters that will be global for all sequence elements."""
    
    def __init__(self, channels=('x', 'y', 'z'), rf_spoil=True, **kwargs):
        """Construct parameter class.

        Parameters
        ----------
        channels : tuple, optional
            Length 3 tuple with 'x', 'y', and 'z'.  All other code should use channels to
            set the readout, phase-encode, and slice-select directions.  Note that this is
            mostly to be used to force change in axis during scanning, if you want to
            allow rotations on the scanner, then just leave it to the default, and the
            intperretor will rotate for you.
            By default ('x', 'y', 'z')
        rf_spoil : bool, optional
            Is rf spoiling being used, by default True
        """
        
        opt_args = {
                'grad_unit': 'mT/m',
                'slew_unit': 'T/m/s',
                'rf_ringdown_time': 30e-6,
                'rf_dead_time': 100e-6,
                'adc_dead_time': 10e-6,
                'B0': 2.89,
                'adc_samples_limit': 8192,
            }

        # Add any Opts arguments from kwargs to the self Opts
        opt_args.update(
            {key: val for key, val in kwargs.items() if key in Opts.default.__dict__},
        )

        self.system = Opts(**opt_args)
        self.channels = channels
        self.rf_spoil = rf_spoil

        self.rf_spoil_idx = 0
        self.rf_spoil_phase = 0


class PSeq_Base:
    """Base class and common functions for pyPulseq helper classes."""

    def __init__(self, pparams):
        """Construct base class.
        """
        self.pparams = pparams

        self.seq = pp.Sequence(system=self.pparams.system)
        self.track_time = 0

    def increment_rf_spoiling(self):
        """
        Increment the rf spoiling phase.

        rf_spoil_phase should then be used any time an RF or ADC even is played.

        TODO
        ----
        * Add support for different style of rf spoiling (mainly random, or list based)

        """
        self.pparams.rf_spoil_idx += 1
        self.pparams.rf_spoil_phase = (117 * np.pi / 180) * self.pparams.rf_spoil_idx * self.pparams.rf_spoil_idx / 2
        self.pparams.rf_spoil_phase = self.pparams.rf_spoil_phase % (2 * np.pi)

    def reinit_seq(self):
        """Reset the sequence."""
        self.pparams.rf_spoil_idx = 0
        self.pparams.rf_spoil_phase = 0
        self.track_time = 0
        self.seq = pp.Sequence(system=self.pparams.system)

    def add_pseq_to_self(self, pseq_in, *args, **kwargs):
        """Add a sequence element to this sequence.

        Basically call add the pseq_in component to this current self.seq

        Parameters
        ----------
        pseq : Pseq_Base derived
            The sequence element we will add to this sequence

        Returns
        -------
        float
            Duration [seconds] of the added component.
        """
        return pseq_in.add_to_seq(self, *args, **kwargs)

    def add_delay(self, delay):
        """Add a simple delay to the sequence.

        Parameters
        ----------
        delay : float
            Time of the selay in [seconds]
        """
        self.seq.add_block(pp.make_delay(delay))

    def get_seq_time(self):
        """Return the time in seconds of the entire sequence."""
        return self.seq.duration()[0]

    def add_dummy_adc(self):
        """
        Add a short ADC to the sequence (1.0msec).

        This is for Skope sequences with no ADC, I used to get an error when there was no
        ADC, but I dont think this happens anymore.  It still might be useful to make sure
        the scanner actually records something during the sequence.
        """
        self.seq.add_block(
            pp.make_adc(
                num_samples=100,
                duration=100 * 4e-6,
                delay=self.pparams.system.adc_dead_time,
                system=self.pparams.system,
            ),
            pp.make_delay(1e-3),
        )

    def build_blocks(self):
        """Abstract base class for assembling sequence elements."""
        raise NotImplementedError

    def make_default_seq(self, *args, **kwargs):
        """Make a dummy sequence of this class."""
        self.reinit_seq()
        all_blocks = self.build_blocks(*args, **kwargs)
        for block in all_blocks:
            self.seq.add_block(*block)

    def add_to_seq(self, pseq0, *args, **kwargs):
        """
        Add this sequence component to a parent sequence `pseq0`.

        See each components build_blocks function for argument list.

        TODO
        ----
        * Consider if the rf_spoil blanket application will be OK.  i.e. is there ever a
          case where some rf or adc elements should not have the same phase?  If so how
          can we handle it?  We can move it back to build_blocks like it used to be.
        * I think I should move these back to individual classes so that the arguments can
          be more easily determined for IDEs, right now people need to go find the right
          build_blocks to learn arguments.  (We can still keep the build_blocks usage)
        """
        all_blocks = self.build_blocks(*args, **kwargs)
        for block in all_blocks:
            # Set rf spoiling if needed to all rf and adc blocks
            if pseq0.pparams.rf_spoil:
                for bb in block:
                    if bb.type in ['rf', 'adc']:
                        bb.phase_offset = pseq0.pparams.rf_spoil_phase
            
            # Add block to sequence
            pseq0.seq.add_block(*block)

        total_dur = self.get_duration(all_blocks=all_blocks)
        pseq0.track_time += total_dur

        return total_dur

    def get_duration(self, all_blocks=None, *args, **kwargs):
        """Get duration of sequence component."""
        if all_blocks is None:
            all_blocks = self.build_blocks(*args, **kwargs)

        dur = 0
        for block in all_blocks:
            dur += pp.calc_duration(*block)

        return dur

    def check(self, report_time = True):
        """Run pypulseq checks on the sequence, and report scan time."""
        if report_time:
            seq_time = self.get_seq_time()
            print('Seq Time:', datetime.timedelta(seconds=round(seq_time)))

        # Error "p_vit = ..." is from calculating trajectory during empty ADC, seems fine
        ok, error_report = self.seq.check_timing()
        if ok:
            print("Timing check passed successfully")
        else:
            print("Timing check failed. Error listing follows:")
            [print(e) for e in error_report]


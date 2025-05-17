Sequence generation files for variouss forms of GIRFs

# Getting Started
A simple thin-slice measurement can be found in `make_simple_phantom_girf.ipynb`, that will demonstrate how to make a 5-slice, chirp-based measurement that can be acquired in a phantom, using multiple chirp waveforms (with slightly different parameters and prephasers).  It will acquire data on each gradient axis with no spatial encoding.


# pypulseq Helpers
To modularize the code, most pypulseq functionality is packaged inside of object seen in the `helpers/` directory.  

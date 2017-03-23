# ScS_reverb_setup

This directory contains all scripts related to setting up simulation parameters.
No data should live here. Might have simple plotting scripts to make sure that
parameters are being set up correctly.


**make_long_migrate_wave.py** Creates an h5 file with wavelets from 
discontinuities in the mantle. It does this by windowing the 670 reverberation
wavelet from a synthetic and shifting it around. It has the option to shift
depth phases in order to simulate deep events in which the 660 wavelet 
interferes with the zeroth-order phase

**long_migrate.py** performs the migration on a data trace. It uses the synthetic
migrations made my 'make_long_migrate_wave.py'

**make\_lookup\_table** makes an HDF5 file with lookup values for ScS reverberation
traveltimes.


**mineos\_perturb\_220** creates a suite of 1D models which are equal to PREM below
the 220km discont. The perturbation shifts S wave velocity above the 220 by any 
factor you choose. The model is written to be used for mineos.

**perturb\_depth** makes mineos source files that have different depths. This
is to explore the effect of source depth on depth phase separation.

**dip\_perturb\_mineos** makes source files with perturbations to dip. This is to 
better match ScSn and sScSn waveform amplitudes.

**make\_wave\_glossary** makes an HDF5 file with first order reverberation waveforms
included.

**invert\_gf** is a simple Gm=d inversion for greens functions.

**reverb\_greens\_functions** is used to make greens functions for arbitrary
mid mantle first order reverberations.

**mask\_zeroth\_order** is used to remove the zeroth order ScSn waveforms
from the data by subtracting the best fit synthetic

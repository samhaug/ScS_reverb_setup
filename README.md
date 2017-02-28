# ScS_reverb_setup

This directory contains all scripts related to setting up simulation parameters.
No data should live here. Might have simple plotting scripts to make sure that
parameters are being set up correctly.


mineos_perturb_220 creates a suite of 1D models which are equal to PREM below
the 220km discont. The perturbation shifts S wave velocity above the 220 by any 
factor you choose. The model is written to be used for mineos.

perturb_depth makes mineos source files that have different depths. This
is to explore the effect of source depth on depth phase separation.

dip_perturb_mineos makes source files with perturbations to dip. This is to 
better match ScSn and sScSn waveform amplitudes.

make_lookup_table makes an HDF5 file with lookup values for ScS reverberation
traveltimes.

make_wave_glossary makes an HDF5 file with first order reverberation waveforms
included.

Soon to be greens function source generator

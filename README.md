# ScS_reverb_setup

This directory contains all scripts related to setting up simulation parameters.
No data should live here. Might have simple plotting scripts to make sure that
parameters are being set up correctly.

**Reverberation workflow**

**Steps 1-7**
Find optimal perturbed velocity model and perturbed source depth for modeling
reverberations.

**Steps 8-10**
Make green's functions for higher order reverberations and migrate.

**1.**
Identify source/network. Remove dirty seismograms. Use tools in /mineos\_bin to
make input files for mineos. Need short.site, short.sitechan, and source file
made from NDK copied from Wilber.

**2.**
Simulate event with PREM background model using mineos.

**3.**
Use tools in /ScS\_reverb\_plot to compare data to synthetics. Use verify\_stack or
verify\_section. Look for time delays in ScS\_n when aligned on ScS\_2.

**4.**
Use vel\_perturb\_mineos to make perturbed PREM input files for mineos. These
models will have perturbed S wave velocity above 220km. Simulate the event again
for each perturbed velocity model.

**5.**
Repeat step 3. The goal is to find the best perturbed velocity model for predicting
ScS\_n traveltimes.

**6.**
Once optimal perturbed velocity model is found look for depth phase separation.
The next step is to use depth\_perturb\_mineos to make source files with 
perturbed depths. Simulate event with each perturbed depth using the best fit
velocity model from step 4.

**7.**
Repeat step 3, optional to use verify\_depth\_phase to find best depth. 
You now should have an optimal perturbed velocity and source depth model for the
source/network in step 1. This model should predict reverberation traveltimes.

**optional 8. (If source depth > 400km**
Use the best fit velocity model and generate an earthquake of source depth 350km.
These reverberations can be identified and shifted.

**8.**
Use make\_lookup\_table in conjunction with best fit synthetic and best fit 
velocity model to make hdf5 file of traveltimes. This is used for making the 
higher order traveltime green's functions. make\_10km\_taup turns a mineos
input velocity file into a taup lookup .tvel with 10km discontinuity spacing.

**9.**
Use make\_long\_migrate\_wave to construct reverberation green's functions. This
function is where most of the work is done. It expects a traveltime lookup table
and a single trace. It will generate an hdf5 file of greens function reverberations.
If optional 8 was done, strip 670 reverberations from 350km depth and shift them
to their expected arrival times.

**10.**
Use long\_migrate to perform the migration. It expects the greens function
library and a data trace. 



**Scripts in this directory**

**make\_long\_migrate\_wave.py** Creates an h5 file with wavelets from 
discontinuities in the mantle. It does this by windowing the 670 reverberation
wavelet from a synthetic and shifting it around. It has the option to shift
depth phases in order to simulate deep events in which the 660 wavelet 
interferes with the zeroth-order phase

**long_migrate.py** performs the migration on a data trace. It uses the synthetic
migrations made my 'make\_long\_migrate\_wave.py'

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

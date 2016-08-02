# OCL-scripts

Scripts by Jørgen E. Midtbø for stuff, such as data analysis, at the Oslo Cyclotron Lab.
Questions or comments are welcome: j.e.midtbo {at] fys.uio.no

## The OCL calibration library
This is an extensive library of Python functions for aligning and calibrating spectra from SiRi and CACTUS. It should also be perfectly possible to apply them to other types of spectra, potentially with minor tweaks. 

There are example scripts showing how to use them, with corresponding spectrum files.

An example workflow could go something like this:
* Sort data with plain gainshifts.
* Do a coarse alignment of SiRi (EdE) spectra, getting relative gains and shifts (GS). This might be easiest using the peakfinder script (see below). Combine relative GS with previous GS. Resort data.
* Look at the aligned spectra in ROOT. Apply TCut to interesting region(s). Resort data with TCuts applied (see user_sort.cpp).
* Improve the SiRi alignment using only the interesting regions. Get relative GS, combine with previous, resort.
* Add aligned spectra together to have good statistics on peaks. Fit two peaks in E and dE spectra. Determine physical values of the peaks, calibrate data. The function calibrate_2points_from_list_summed_backs() may be helpful. Make relative GS file, combine with previous, resort.
* Start with the CACTUS (NaI) spectra. Align all spectra. Get relative GS, combine with previous, resort.
* Sum spectra to get nice peaks. Fit peaks. Determine physical values, calculate gain and shift. Combine with previous, resort. 
* Do time gating and other coincidence stuff.
* Proceed with $E_x-E_\gamma$ matrix.


## The peakfinder script
peaks2D_JEMmod.C is a script by Alexander Bürger that is very nice for determining the locations of peaks in 2D SiRi spectra. I have modified it slightly such that pressing "u" moves you to the next spectrum in line, while pressing "space" gives you another shot at the same spectrum. 

After obtaining two .csv files with coordinates for two different peaks, for each of the 64 front and back detectors, the function calibrate_2points() from the OCL calibration library may be used to make a gainshifts file. The .csv files have to have only one entry per detector, so duplicates should be removed manually first. Also, calibrate_2points() assumes that you know which physical energies the peaks correspond to, so you have to make a choice (although it can be arbitrary if you're only aligning). This requires a set of Qkinz files named in a specific way, see example folder in this repository. Qkinz can be downloaded from https://github.com/oslocyclotronlab/Qkinz .
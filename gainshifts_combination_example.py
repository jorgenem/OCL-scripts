# Made by J{\o}rgen E. Midtb{\o}, University of Oslo, 20160406
# j.e.midtbo@fys.uio.no

from OCL_calibration_library import *
# Example script for combining two gainshift files into a new, composite file. 
# Requires the example gainshift files in the git repository. (Note the formatting
# of the RELATIVE files.)

# To update EdE gains and shifts
make_composite_gainshifts_ede("gainshifts-0-plain.dat", "gainshifts-RELATIVE-EdE.dat", "gainshifts-1-EdE_composite.dat")

# To update NaI gains and shifts:
make_composite_gainshifts_nai("gainshifts-0-plain.dat", "gainshifts-RELATIVE-NaI.dat", "gainshifts-1-NaI_composite.dat")
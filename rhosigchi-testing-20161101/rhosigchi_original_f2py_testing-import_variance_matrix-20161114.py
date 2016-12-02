import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/OCL-scripts')
import pyma_v1 as pyma
sys.path.insert(0,'/home/jorgenem/gitrepos/oslo-method-software/prog')
import rhosigchi_f2py_importvar
import numpy as np 

print rhosigchi_f2py_importvar.__doc__

fg, calib, Eg_range, Ex_range = pyma.read_mama('firstgen-jem-20161011-Nexbins196-Nstat500.m')
fgv, tmp1, tmp2, tmp3 = pyma.read_mama('fgvar-jem-20161011-Nexbins196-Nstat500.m')
# print calib
# print Ex_range, Eg_range
import matplotlib.pyplot as plt

Eg_min = 1000
Ex_min = 4000
Ex_max = 8000


# # Rebinning matrices to get right format, should be 120 keV bins
# N = 65
# fg = fg[0:193,:]
# fgv = fgv[0:193,:]
# Ex_range = Ex_range[0:193]
# Eg_range = Eg_range[0:193]
# fg, Ex_range_rebinned = pyma.rebin_and_shift(fg, Ex_range, N, rebin_axis=0)
# fg, Eg_range_rebinned = pyma.rebin_and_shift(fg, Eg_range, N, rebin_axis=1)
# fgv, tmp = pyma.rebin_and_shift(fgv, Ex_range, N, rebin_axis=0)
# fgv, tmp = pyma.rebin_and_shift(fgv, Eg_range, N, rebin_axis=1)
# calib = [Eg_range_rebinned[0], Eg_range_rebinned[1]-Eg_range_rebinned[0], 
# 		 Ex_range_rebinned[0], Ex_range_rebinned[1]-Ex_range_rebinned[0]]
# print calib


# 20161114: Test of just sending in original matrix and letting rhosigchi do rebinning with ELASTIC()


# plt.matshow(fg, origin='lower')
# plt.show()
# sys.exit(0)

# sys.exit(0)

# Reformat: Put fg and fgv into 512x512 matrices to comply with rhosigchi:
fg_reformatted = np.zeros((512,512))
fg_reformatted[0:fg.shape[0],0:fg.shape[1]] = fg
fgv_reformatted = np.zeros((512,512))
fgv_reformatted[0:fgv.shape[0],0:fgv.shape[1]] = fgv
fg_reformatted = fg_reformatted.T 	# Fortran uses 
fgv_reformatted = fgv_reformatted.T # column-major order






# sys.exit(0)
rho, T = rhosigchi_f2py_importvar.rhosigchi(fg_reformatted, fgv_reformatted, calib, Eg_min, Ex_min, Ex_max)

plt.figure()
plt.plot(rho,label='rho new')
plt.plot(T,label='T new')

import rhosigchi_f2py_origvar_pythonoutput
rho_orig, T_orig = rhosigchi_f2py_origvar_pythonoutput.rhosigchi(fg_reformatted, calib, Eg_min, Ex_min, Ex_max)

plt.plot(rho_orig,label='rho orig')
plt.plot(T_orig,label='T orig')

plt.yscale('log')
plt.legend()
plt.show()

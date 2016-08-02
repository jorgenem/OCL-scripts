from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from matplotlib.colors import LogNorm

filenames = [
'alfna-attempted_removing_1p2MeV_Al-20160518.m',
'alfna-unfolded-20160518.m',
'alfna-20160518.m',
'alfna-1stgen-20160518.m',
'alfna-attempted_removing_2p2MeV_Al-20160518.m'
]

# filename = 'alfna-attempted_removing_2p2MeV_Al-20160518.m'
for filename in filenames:
	matrix = np.genfromtxt(filename, skip_header=10, skip_footer=3)
	datafile = open(filename, 'r')
	calibration_line = datafile.readlines()[6].split()
	a = [float(calibration_line[2][:-1])*1e-3, float(calibration_line[3][:-1])*1e-3, float(calibration_line[5][:-1])*1e-3, float(calibration_line[6][:-1])*1e-3]
	Nx = len(matrix[0,:])
	Ny = len(matrix[:,0])
	x_array = np.linspace(a[0], a[0]+a[1]*(Nx+1), Nx+1)
	y_array = np.linspace(a[2], a[2]+a[3]*(Ny+1), Ny+1)
	
	plt.figure()
	plt.pcolormesh(x_array, y_array, matrix, norm=LogNorm(vmin=0.01, vmax=matrix.max()), cmap='gist_rainbow_r')
	plt.colorbar()
	plt.xlabel('$E_\gamma$ [MeV]', fontsize=14)
	plt.ylabel('$E_x [MeV]$', fontsize=14)
	plt.ylim([0,12])
	plt.xlim([0,12])
	
	plt.savefig(filename[:-2]+'_pyplot.png')
	# plt.show()







	datafile.close()
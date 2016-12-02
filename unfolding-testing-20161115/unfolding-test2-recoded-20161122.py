from __future__ import division

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/OCL-scripts')
import pyma_v1 as pyma

# Import raw mama matrix
raw, calib, Eg_range, Ex_range = pyma.read_mama('alfna-20160518.m')
R, tmp1, tmp2, tmp3 = pyma.read_mama('response-20161122.m')

# raw_reformat = np.zeros((4096,2048))
# raw_reformat[0:raw.shape[0],0:raw.shape[1]] = raw
# # raw_reformat[0:raw.T.shape[0],0:raw.T.shape[1]] = raw.T # Test with transposition
# # raw_reformat[10,100] = 100
# print raw_reformat.max()

# plt.matshow(raw_reformat,origin='lower')
# plt.show()

# plt.matshow(R,origin='lower')
# plt.show()

# sys.exit(0)

# Run unfolding
# rawmat = np.zeros((4096,2048))
calib = calib[0:2]
print "calib = ",calib

# Set limits for excitation and gamma energy bins to be considered for unfolding
IEx_low = 0
IEx_high = 250
# IEx_high = raw.shape[0]
IEg_low = 0
IEg_high = 1000
# IEg_high = raw.shape[1]
Nit = 35

# Make masking array to cut away noise below Eg=Ex+dEg diagonal
# Define cut   x1    y1    x2    y2
cut_points = [ 72,   5,  1050,  257]
def line(x, points):
    a = (points[3]-points[1])/float(points[2]-points[0])
    b = points[1] - a*points[0]
    print "a = {}, b = {}".format(a,b)
    return a*x + b
i_array = np.linspace(0,len(Ex_range)-2,len(Ex_range)-1).astype(int) # Ex axis 
j_array = np.linspace(0,len(Eg_range)-2,len(Eg_range)-1).astype(int) # Eg axis
i_mesh, j_mesh = np.meshgrid(i_array, j_array, indexing='ij')
mask = np.where(i_mesh > line(j_mesh, cut_points), 1, 0)



rawmat = (raw*mask)[IEx_low:IEx_high,IEg_low:IEg_high].T # Has to have axis ordering (Ex,Eg) to fit f = R*u matrix product
Ndof = mask[IEx_low:IEx_high,IEg_low:IEg_high].sum()

unfoldmats = np.zeros((rawmat.shape[0],rawmat.shape[1],Nit))
foldmats = np.zeros((rawmat.shape[0],rawmat.shape[1],Nit))
chisquares = np.zeros(Nit)
R = R[IEg_low:IEg_high,IEg_low:IEg_high]

# plt.plot(R[300,:])
# plt.show()

# sys.exit(0)

unfoldmats[:,:,0] = rawmat # First approximation to unfolded spectrum is raw(=folded) spectrum
for iteration in range(Nit-1):
    if iteration>0:
        unfoldmats[:,:,iteration] = unfoldmats[:,:,iteration-1]+(unfoldmats[:,:,0]-foldmats[:,:,iteration-1])
    foldmats[:,:,iteration] = np.dot(R,unfoldmats[:,:,iteration])
    chisquares[iteration] = pyma.div0(np.power(foldmats[:,:,iteration]-rawmat,2),np.where(rawmat>4,rawmat,4)).sum() / Ndof
    print "Iteration = {}, chisquare = {}".format(iteration,chisquares[iteration])
    if iteration < 3 or iteration > 25:
        unfoldmat = np.copy(unfoldmats[:,:,iteration].T)
        # unfoldmat[unfoldmat<=0] = NaN
        # plt.subplot(2,2,1)
        # plt.pcolormesh(unfoldmat,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=0, vmax=270))
        # plt.subplot(2,2,2)
        # plt.pcolormesh(foldmats[:,:,iteration].T,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=0, vmax=270))
        # plt.subplot(2,2,3)
        # plt.pcolormesh(np.power(rawmat-foldmats[:,:,iteration],2).T,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=0, vmax=270))
        # plt.colorbar()
        # plt.subplot(2,2,4)
        # plt.pcolormesh(rawmat.T,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=0, vmax=270))
        # plt.show()

pyma.write_mama(unfoldmats[:,:,15].T, 'pyma_unfolded_test-20161129.m', Eg_range, Ex_range[IEx_low:IEx_high])
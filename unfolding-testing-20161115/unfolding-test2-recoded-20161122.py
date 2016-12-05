from __future__ import division

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# from scipy.ndimage.filters import gaussian_filter, gaussian_filter1d
from scipy.stats import norm

import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/OCL-scripts')
import pyma_v1 as pyma


def i_from_E(E, E_range):
    # Function which returns the index of the E_range value closest to given E
    where_array = np.where(E_range > E)[0]
    # print where_array, len(where_array)
    if len(where_array) > 0:
        i = where_array[0]
        if np.abs(E_range[i]-E) > np.abs(E_range[i-1]-E):
            i -= 1
    else:
        i = len(E_range)-1
    return i

def shift_and_smooth(array, E_range, FWHM, shift):
    # Takes an array of counts, shifts it (downward only!) with energy 'shift'
    # and smooths it with a gaussian of specified 'FWHM'
    i_shift = i_from_E(shift, E_range)
    indices = (np.linspace(0,len(array)-1,len(array)).astype(int)+i_shift)
    indices = np.trim_zeros(np.where(indices < len(array), indices, 0), trim='b')
    # print indices
    array_out = array[indices]
    if len(array_out) < len(array):
        array_out = np.append(array_out, np.zeros(len(array)-len(array_out)))

    # Apply smoothing
    array_out = gaussian_filter(array_out,FWHM/2.355)
    return array_out

def shift_and_smooth2D(array, E_range, FWHM, p, shift, smoothing=True):
    # Takes an array of counts, shifts it (downward only!) with energy 'shift'
    # and smooths it with a gaussian of specified 'FWHM'.
    # This version is vectorized to shift, smooth and scale all points
    # of 'array' individually, and then sum together and return.

    # The arrays from resp.dat are missing the first channel.
    p = np.append(0, p) 
    FWHM = np.append(0, FWHM)

    a1 = (E_range[1]-E_range[0]) # bin width
    N = len(array)
    # Shift is the same for all energies 
    if shift == "annihilation":
        # For the annihilation peak, all channels should be mapped on E = 511 keV. Of course, gamma channels below 511 keV,
        # and even well above that, cannot produce annihilation counts, but this is taken into account by the fact that p
        # is zero for these channels. Thus, we set i_shift=0 and make a special indices_shifted array to map all channels of
        # original array to i(511). 
        i_shift = 0 
    else:
        i_shift = i_from_E(shift, E_range) - i_from_E(0, E_range) # The number of indices to shift by
    print i_shift
    # indices = (np.linspace(0,len(array)-1,len(array)).astype(int)+i_shift) # Index array for shifted array
    # indices = np.trim_zeros(np.where(indices < len(array), indices, 0), trim='b') # Make sure nothing wraps around and comes down from the top
    indices_original = np.linspace(i_shift, len(array)-1, len(array)-i_shift).astype(int) # Index array for original array, truncated to shifted array length
    if shift == "annihilation": # If this is the annihilation peak then all counts should end up with their centroid at E = 511 keV
        indices_shifted = (np.ones(N-i_from_E(511, E_range))*i_from_E(511, E_range)).astype(int)
    else:
        indices_shifted = np.linspace(0,len(array)-i_shift-1,len(array)-i_shift).astype(int) # Index array for shifted array

    # Transform energy array from edge to middle-bin
    E_range_middlebin = (E_range[0:N]+E_range[1:N+1])/2

    if smoothing:
        matrix = np.zeros((N,N))
        indices_original = indices_original[np.where(FWHM[indices_original]>0)] # Filter out channels with zero width, we don't want those
        for i in indices_original: # i is the energy channel in unshifted array
            try:
                matrix[i] = array[i]*p[i]*norm.pdf(E_range_middlebin, loc=Eg_range_middlebin[indices_shifted[i]], scale=FWHM[i]/2.355) \
                * a1 # Multiplying by bin width to preserve number of counts
                            # TODO: Figure out how FWHM relates to bin width
            except IndexError:
                pass

        # plt.matshow(matrix)
        # plt.show()
        array_out = matrix.sum(axis=0)
    else:
        array_out = np.zeros(N)
        # for i in range(N):
        #     try:
        #         array_out[i-i_shift] = array[i] #* p[i+1]
        #     except IndexError:
        #         pass

        # Instead of above, vectorizing:
        array_out[indices_shifted] = p[indices_original]*array[indices_original]

    return array_out


# Import raw mama matrix
raw, calib, Eg_range, Ex_range = pyma.read_mama('alfna-filled_negative-20161205.m')
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




# === Run iterative unfolding ===



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

# pyma.write_mama(unfoldmats[:,:,15].T, 'pyma_unfolded_test-fn-20161205.m', Eg_range, Ex_range[IEx_low:IEx_high])



# === Start compton subtraction method ===

unfolded = unfoldmats[:,:,25].T
N = unfolded.shape[1]
print "N = {}".format(N)

# raw, calib_raw, Eg_raw, Ex_raw = pyma.read_mama('alfna-20160518.m')
# raw, calib_raw, Eg_raw, Ex_raw = pyma.read_mama('alfna-filled_negative-20161205.m')
raw = raw[0:250, 0:1000] # TODO: Should probably mask this to remove noise
print "raw.shape = {}, unfolded.shape = {}".format(raw.shape, unfolded.shape)


# Read resp.dat
resp = np.loadtxt('resp.dat',skiprows=4)

# TODO: Write gauss smoothing function. Take in FWHM, (position of maximum?). Remember to normalize. Compare to Magne's code. Vectorize to do all Ex at once?

i_Eg = 200

# Here is one of the actual unfolded spectra for testing after verifying it works on a simple peak:
orig = unfolded[i_Eg,:]
raw_array = raw[i_Eg,:]

plt.figure(1)
plt.subplot(2,1,1)
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, orig, label='Iterated')
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, raw_array, label='Raw')
plt.legend()
# orig = np.zeros(N)
# orig[i_from_E(5000,Eg_range)] = 400
# orig[i_from_E(7000,Eg_range)] = 300

FWHM = resp[:,1]
eff = resp[:,2]
pf = resp[:,3]
pc = resp[:,4]
ps = resp[:,5]
pd = resp[:,6]
pa = resp[:,7]
Eg_range_middlebin = (Eg_range[0:N]+Eg_range[1:N+1])/2
# fe = orig * np.append(0, pf[0:N-1])
fe = shift_and_smooth2D(orig, Eg_range, 0.5*FWHM, pf, shift=0, smoothing=True)
se = shift_and_smooth2D(orig, Eg_range, 0.5*FWHM, ps, shift=511, smoothing=True)
de = shift_and_smooth2D(orig, Eg_range, 0.5*FWHM, pd, shift=1022, smoothing=True)
a  = shift_and_smooth2D(orig, Eg_range, FWHM, pa, shift="annihilation", smoothing=True)

w = se + de + a 

v = fe + w
c = raw_array - v

c_s = shift_and_smooth2D(c, Eg_range, FWHM, np.ones(len(FWHM)), shift=0, smoothing=True)

u = pyma.div0((raw_array - c_s - w), np.append(0,pf)[0:N])
U = pyma.div0(u, np.append(0, eff)[0:N])

# compton_check = orig * np.append(0, pc[0:N-1]) # All the counts that go into compton background, but at their original channels, just to check that we're not losing any counts



plt.subplot(2,1,2)
plt.plot(Eg_range_middlebin, v, label='v (=p*f+se+de+a)')
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, raw_array, label='Raw')
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, c, label='Compton (r-v)')
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, c_s, label='Compton smoothed')
plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, u, label='Unfolded')
# plt.step((Eg_range[0:N]+Eg_range[1:N+1])/2, U, label='Unfolded, eff corrected')
plt.legend()
plt.show()

print "Sums: orig = {}, v = {}, v+compton = {}, v+smoothed compton = {}".format(orig.sum(), v.sum(), (v+c).sum(), (v+c_s).sum())
sys.exit(0)
# both = pf*orig + ps*smoothed

# print 





plt.figure()
plt.step(Eg_range_middlebin, orig,label='orig')
plt.step(Eg_range_middlebin, smoothed,label='modified')
plt.legend()
plt.show()
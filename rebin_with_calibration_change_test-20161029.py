from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 

def rebin_and_shift(array, E_range, N_final, rebin_axis=0):
	# Function to rebin an M-dimensional array either to larger or smaller binsize.
	# Written by J{\o}rgen E. Midtb{\o}, University of Oslo, j.e.midtbo@fys.uio.no, github.com/jorgenem
	# Latest change made 20161029.

	# Rebinning is done with simple proportionality. E.g. for down-scaling rebinning (N_final < N_initial): 
	# if a bin in the original spacing ends up between two bins in the reduced spacing, 
	# then the counts of that bin are split proportionally between adjacent bins in the 
	# rebinned array. 
	# Upward binning (N_final > N_initial) is done in the same way, dividing the content of bins
	# equally among adjacent bins.

	# Technically it's done by repeating each element of array N_final times and dividing by N_final to 
	# preserve total number of counts, then reshaping the array from M dimensions to M+1 before flattening 
	# along the new dimension of length N_initial, resulting in an array of the desired dimensionality.

	# This version (called rebin_and_shift rather than just rebin) takes in also the energy range array (lower bin edge)
	# corresponding to the counts array, in order to be able to change the calibration. What it does is transform the
	# coordinates such that the starting value of the rebinned axis is zero energy. This is done by shifting all
	# bins, so we are discarding some of the eventual counts in the highest energy bins. However, there is usually a margin.
	
	N_initial = array.shape[rebin_axis] # Initial number of counts along rebin axis

	# Repeat each bin of array Nfinal times and scale to preserve counts
	array_rebinned = array.repeat(N_final, axis=rebin_axis)/float(N_final)

	if E_range[0] < 0 or E_range[1] < E_range[0]:
		print "Erorr in function rebin_and_shift(): Negative zero energy is not supported. (But it should be relatively easy to implement.)"
		sys.exit(0)

	# Calculate number of extra slices in Nf*Ni sized array required to get down to zero energy
	n_extra = int(np.ceil(N_final * (E_range[0]/(E_range[1]-E_range[0]))))
	# Append this matrix of zero counts in front of the array
	indices_append = np.array(array_rebinned.shape)
	indices_append[rebin_axis] = n_extra
	array_rebinned = np.append(np.zeros(indices_append), array_rebinned, axis=rebin_axis)
	array_rebinned = np.split(array_rebinned, [0, N_initial*N_final], axis=rebin_axis)[1]
	indices = np.insert(array.shape, rebin_axis, N_final) # Indices to reshape to
	array_rebinned = array_rebinned.reshape(indices).sum(axis=(rebin_axis+1)) 
	E_range_shifted_and_scaled = np.linspace(0, E_range[-1]-E_range[0], N_final)
	return array_rebinned, E_range_shifted_and_scaled

np.random.seed(2)
N_initial = 26
E_range_initial = np.linspace(15,415,N_initial)
counts_initial = np.random.uniform(size=N_initial)

plt.figure(0)
plt.subplot(1,2,1)
plt.step(E_range_initial, counts_initial)

N_final = 60
counts_final, E_range_final = rebin_and_shift(counts_initial, E_range_initial, N_final)

plt.subplot(1,2,2)

plt.step(E_range_final, counts_final)

print "E_range_initial =",E_range_initial
print "E_range_final =",E_range_final


plt.show()
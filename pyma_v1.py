from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from matplotlib.colors import LogNorm


def read_mama(filename):
	# Reads a MAMA matrix file and returns the matrix as a numpy array, 
	# as well as a list containing the four calibration coefficients
	# (ordered as [bx, ax, by, ay] where Ei = ai*channel_i + bi)
	# and 1-D arrays of calibrated x and y values for plotting and similar.
	matrix = np.genfromtxt(filename, skip_header=10, skip_footer=1)
	datafile = open(filename, 'r')
	calibration_line = datafile.readlines()[6].split()
	a = [float(calibration_line[2][:-1]), float(calibration_line[3][:-1]), float(calibration_line[5][:-1]), float(calibration_line[6][:-1])]
	# TODO: INSERT CORRECTION FROM CENTER-BIN TO LOWER EDGE CALIBRATION HERE.
	# MAKE SURE TO CHECK rebin_and_shift() WHICH MIGHT NOT LIKE NEGATIVE SHIFT COEFF.
	# (alternatively consider using center-bin throughout, but then need to correct when plotting.)
	Nx = len(matrix[0,:])
	Ny = len(matrix[:,0])
	# print Nx, Ny
	x_array = np.linspace(a[0], a[0]+a[1]*(Nx), Nx) # BIG TODO: This is probably center-bin calibration, 
	y_array = np.linspace(a[2], a[2]+a[3]*(Ny), Ny) # and should be shifted down by half a bin?
																									# Update 20171024: Started changing everything to lower bin edge,
																									# but started to hesitate. For now I'm inclined to keep it as
																									# center-bin everywhere. 
	return matrix, a, x_array, y_array # TODO: Change this so it returns the y array first - this is axis 0

def write_mama(matrix, filename, Egamma_range, Ex_range):
	import time
	outfile = open(filename, 'w')

	# Write mandatory header:
	# outfile.write('!FILE=Disk \n')
	# outfile.write('!KIND=Spectrum \n')
	# outfile.write('!LABORATORY=Oslo Cyclotron Laboratory (OCL) \n')
	# outfile.write('!EXPERIMENT=pyma \n')
	# outfile.write('!COMMENT=none|RE:alfna-20FN:RN:UN:FN:RN: \n')
	# outfile.write('!TIME=DATE:'+time.strftime("%d-%b-%y %H:%M:%S", time.localtime())+'   \n')
	# outfile.write('!CALIBRATION EkeV=6, %12.6E, %12.6E, 0.000000E+00, %12.6E, %12.6E, 0.000000E+00 \n' %(Egamma_range[0], (Egamma_range[1]-Egamma_range[0]), Ex_range[0], (Ex_range[1]-Ex_range[0])))
	# outfile.write('!PRECISION=16 \n')
	# outfile.write('!DIMENSION=2,0:%4d,0:%4d \n' %(len(matrix[:,0]), len(matrix[0,:])))
	# outfile.write('!CHANNEL=(0:%4d,0:%4d) \n' %(len(matrix[:,0]), len(matrix[0,:])))
	header_string ='!FILE=Disk \n'
	header_string +='!KIND=Spectrum \n'
	header_string +='!LABORATORY=Oslo Cyclotron Laboratory (OCL) \n'
	header_string +='!EXPERIMENT=pyma \n'
	header_string +='!COMMENT=none|RE:alfna-20FN:RN:UN:FN:RN: \n'
	header_string +='!TIME=DATE:'+time.strftime("%d-%b-%y %H:%M:%S", time.localtime())+'   \n'
	header_string +='!CALIBRATION EkeV=6, %12.6E, %12.6E, 0.000000E+00, %12.6E, %12.6E, 0.000000E+00 \n' %(Egamma_range[0], (Egamma_range[1]-Egamma_range[0]), Ex_range[0], (Ex_range[1]-Ex_range[0]))
	header_string +='!PRECISION=16 \n'
	header_string +='!DIMENSION=2,0:%4d,0:%4d \n' %(len(matrix[0,:])-1, len(matrix[:,0])-1)
	header_string +='!CHANNEL=(0:%4d,0:%4d) ' %(len(matrix[0,:])-1, len(matrix[:,0])-1)

	footer_string = "!IDEND=\n"

	# Write matrix:
	# matrix.tofile(filename, sep='       ', format="{:14.8E}")
	# matrix.tofile(filename, sep=' ', format="%-17.8E")
	np.savetxt(filename, matrix, fmt="%-17.8E", delimiter=" ", newline="\n", header=header_string, footer=footer_string, comments="")

	outfile.close()



# Remove contaminations:
def remove_contamination_peak(matrix, limits):
	x1, x2, y1, y2 = limits
	Ny, Nx = y2-y1, x2-x1
	Ny_tot, Nx_tot = matrix.shape
	# Subtract linear background for each Ex bin:
	shape_individual = np.zeros((Ny, Nx))
	for iy in range(Ny):
		shape_individual[iy,:] = matrix[iy+y1,x1:x2] - (np.linspace(0,1,Nx)*(matrix[iy+y1,x2]-matrix[iy+y1,x1])/(x2-x1)+matrix[iy+y1,x1])
		# plt.pcolormesh(shape_individual)
	# plt.colorbar()
	# plt.show()

	# Calculate area of shape in each Ex bin:
	area = np.abs(np.sum(shape_individual, axis=1))
	
	# Take the mean shape and normalize to unity:
	shape_avg = np.mean(shape_individual,axis=0)
	shape_avg /= np.sum(shape_avg)
	# plt.plot(shape_avg)
	# plt.show()

	# Subtract the average shape multiplied by the area from each individual shape
	matrix_subtracted = np.copy(matrix)
	for iy in range(Ny):
		# To get the area of the contamination peak in the current iy, we subtract a linear background again and sum the rest.
		# The linear interpolation used to subtract background in this case is based on the mean of the endpoints 
		# at iy-1, iy and iy+1. This is the same way as in MAMA.
		z1 = np.mean(matrix[iy-1+y1:iy+1+y1,x1])
		z2 = np.mean(matrix[iy-1+y1:iy+1+y1,x2])
		z_net = matrix[iy+y1,x1:x2] - (z1 + (z2-z1)/(x2-x1)*np.linspace(0,1,Nx))
		print z_net
		area_current = np.sum(z_net)
		print area_current
		matrix_subtracted[iy+y1,x1:x2] = matrix[iy+y1,x1:x2] - shape_avg*area_current

	plt.figure(10)
	plt.pcolormesh(matrix_subtracted, norm=LogNorm(vmin=0.1, vmax=matrix.max()), cmap='gist_rainbow_r')
	plt.colorbar()
	plt.show()

	return matrix_subtracted

#remove_contamination_peak(matrix, [245, 270, 198, 218])


def slide(Ex, Thres1, Thres2, ThresRatio):
	slide = ThresRatio*Ex
	if slide < Thres1:
		slide = Thres1
	if slide > Thres2:
		slide = Thres2
	return slide


#first_generation_spectrum(matrix, a, N_iterations=20)

def unfold(matrix, a, x_array, y_array, resp_filename):
	return 0


def first_generation_spectrum_test2(matrix, Egamma_range, Ex_range_mat, N_Exbins, Ex_max, dE_gamma, N_iterations=1):
	Ny = len(matrix[:,0])
	Nx = len(matrix[0,:])
	# Extract / calculate calibration coefficients
	bx = Egamma_range[0]
	ax = Egamma_range[1] - Egamma_range[0]
	by = Ex_range_mat[0]
	ay = Ex_range_mat[1] - Ex_range_mat[0]

	statistical_or_total = 1
	ThresSta = 430.0
	# AreaCorr = 1
	ThresTot = 	200.000000
	ThresRatio = 0.3000000
	# ExH = 7520.00000
	ExEntry0s = 300.000000
	ExEntry0t = 0.00000000
	apply_area_correction = True


	# Ex_max = 7500 # keV - maximum excitation energy
	# Ex_min = 300 # keV - minimal excitation energy, effectively moving the ground-state energy up because we cannot resolve the low-energy yrast gamma lines. This is weighed up by also using an effective multiplicity which is lower than the real one, again not considering the low-energy yrast gammas.
	# dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
	# Ex_binsize = 40 # keV - bin size that we want on y axis
	# N_Exbins = 120 # Number of excitation energy bins (NB! It will only rebin in whole multiples, so a slight change in N_Exbins might only result in getting some more empty bins on top.)
	# N_Exbins_original = (Ex_max+dE_gamma)/ay # The number of bins between 0 and Ex_max + dE_gamma in the original matrix
	# grouping = int(np.ceil(len(y_array[np.logical_and(0 < y_array, y_array < Ex_max + dE_gamma)])/N_Exbins)) # The integer number of bins that need to be grouped to have approximately N_Exbins bins between Ex_min and Ex_max after compression (rounded up)

	# Make arrays of Ex and Egamma axis values
	Ex_range = np.linspace(by, Ex_max + dE_gamma, N_Exbins)
	Egamma_range = np.linspace(0,Nx-1,Nx)*ax + bx # Range of Egamma values
	
	# Compress matrix along Ex
	#matrix_ex_compressed = matrix[0:int(N_Exbins*grouping),:].reshape(N_Exbins, grouping, Nx).sum(axis=1)
	matrix_ex_compressed = rebin(matrix[0:int((Ex_max+dE_gamma)/Ex_range_mat.max()*Ny),:], N_Exbins, rebin_axis = 0) # This seems crazy. Does it cut away anything at all?
	# print Ny, N_Exbins, N_Exbins_original	
	# plt.pcolormesh(Egamma_range, Ex_range, matrix_ex_compressed, norm=LogNorm(vmin=0.001, vmax=matrix_ex_compressed.max()))
	# plt.matshow(matrix_ex_compressed)
	# plt.colorbar()
	# plt.show()

	# Remove counts in matrix for Ex higher than Ex_max:
	matrix_ex_compressed[Ex_range>Ex_max, :] = 0
	# plt.matshow(matrix_ex_compressed)
	# plt.colorbar()
	# plt.show()


	# ==== Calculate multiplicities: ====

	# Setup meshgrids for making boolean indexing arrays
	Egamma_mesh, Ex_mesh = np.meshgrid(Egamma_range, Ex_range)
	Egamma_max = Ex_range + dE_gamma # Maximal Egamma value for each Ex bin
	Egamma_max_grid = np.meshgrid(np.ones(Nx), Egamma_max)[1]
	if statistical_or_total == 1:
		# Statistical multiplicity calculation (i.e. trying to use statistical/continuum region only)
		slide = np.minimum( np.maximum(ThresRatio*Ex_mesh, ThresTot), ThresSta ) # The sliding lower limit for Egamma integral - sliding between ThresTot and ThresSta.
		# plt.figure(5)
		# plt.plot(slide[:,0])
		# plt.show()
		# sys.exit(0)
		# good_indices = np.where(np.logical_and(slide < Egamma_mesh, Egamma_mesh < Egamma_max_grid) , True, False)
		matrix_ex_compressed_cut = np.where(np.logical_and(slide < Egamma_mesh, Egamma_mesh < Egamma_max_grid), matrix_ex_compressed, 0)
	elif statistical_or_total == 2:
		# Total multiplicity calculation
		# good_indices = np.where(Egamma_mesh < Egamma_max_grid, True, False)
		matrix_ex_compressed_cut = np.where(Egamma_mesh < Egamma_max_grid, matrix_ex_compressed, 0)
	# for i in range(len(good_indices[:,0])):
	# 	print len(good_indices[i,good_indices[i,:]]) # OK, it actually works.
	
	# Cut away counts higher than Egamma = Ex + dE_gamma
	# matrix_ex_compressed_cut = np.where(good_indices, matrix_ex_compressed, 0)
	# plt.figure(1)
	# plt.pcolormesh(Egamma_range, Ex_range, matrix_ex_compressed_cut, norm=LogNorm(vmin=0.01, vmax=matrix_ex_compressed.max()))
	# plt.show()
	# sys.exit(0)

	# Calculate average multiplicity for each Ex channel
	area_matrix_ex_compressed_cut = np.sum(matrix_ex_compressed_cut, axis=1)
	Egamma_average = div0( np.sum(Egamma_mesh * matrix_ex_compressed_cut, axis =1) , area_matrix_ex_compressed_cut )
	if statistical_or_total == 1:
		# Statistical multiplicity - use the effective Ex0 value
		multiplicity = div0( Ex_range - np.maximum( np.minimum(Ex_range - 200, ExEntry0s), 0), Egamma_average)
	elif statistical_or_total == 2:
		# Total multiplicity - use actual Ex0 = 0
		multiplicity = div0( Ex_range, Egamma_average )

	# plt.figure(2)
	# plt.step(Ex_range, multiplicity) # This rises like a straight line from 0 to about 3-4 - seems very right!
	# plt.show()
	# sys.exit(0)

	# Set up dummy first-generation matrix to start iterations, made of normalized boxes:
	H = np.zeros((N_Exbins, Nx))
	for i in range(N_Exbins):
		Ni = len(Egamma_range[Egamma_range<Ex_range[i] + dE_gamma])
		# print Ni
		H[i, Egamma_range < Ex_range[i] + dE_gamma] = 1/Ni
	# print np.sum(H0, axis=1) # Seems to work!

	# Set up normalization matrix N
	area = np.sum(matrix_ex_compressed_cut, axis=1) # Get total number of counts in each Ex bin
	# plt.plot(Ex_range, area)
	# plt.show()

	area_grid = np.tile(area, (N_Exbins, 1)) # Copy the array N_Exbins times down to make a square matrix
	# print area_grid.shape
	multiplicity_grid = np.tile(multiplicity, (N_Exbins, 1)) 
	# print multiplicity_grid.shape
	normalization_matrix = div0(( np.transpose(multiplicity_grid) * area_grid ) , (multiplicity_grid * np.transpose(area_grid) )).T # The transpose gives the right result. Haven't twisted my head around exactly why.
	# normalization_matrix_check = np.zeros((N_Exbins, N_Exbins))
	# for i in range(N_Exbins):
	# 	for j in range(N_Exbins):
	# 		normalization_matrix_check[i, j] = multiplicity[i]*area[j]/(multiplicity[j]*area[i])
	normalization_matrix[np.isnan(normalization_matrix)] = 0
	# plt.matshow(normalization_matrix, origin='lower', norm=LogNorm(vmin=0.01, vmax=normalization_matrix.max()))
	# plt.show()
	# plt.matshow(normalization_matrix_check, origin='lower') # There is a difference of a transposition, check which is the right one
	# plt.show()

	# Set up compression parameters for Egamma axis to be used by H below:
	i_Egamma_max = np.where(Egamma_range > Ex_max+ dE_gamma)[0][0] # Get the maximal allowed gamma energy (need to make H square, thus Egamma <= Ex + dE_gamma, since that's the maximal Ex channel in the compressed matrix)
	# print i_Egamma_max, Egamma_range[i_Egamma_max], N_Exbins, int(i_Egamma_max/N_Exbins)
	# i_Egamma_max = i_Egamma_max + N_Exbins - i_Egamma_max%N_Exbins # Make sure the number of indices is a whole multiple of N_Exbins (rounded up)
	# print i_Egamma_max
	grouping_Egamma = int(np.ceil(i_Egamma_max/N_Exbins))
	# print grouping_Egamma
	# Egamma_range_compressed = Egamma_range[0:i_Egamma_max]*grouping_Egamma
	Egamma_range_compressed = Ex_range

	# plt.matshow(H[:,0:i_Egamma_max])
	# plt.show()
	# H_extended = np.insert(H[:,0:i_Egamma_max], np.linspace(0,i_Egamma_max, N_Exbins - i_Egamma_max%N_Exbins), H[:,(np.linspace(0,i_Egamma_max, N_Exbins - i_Egamma_max%N_Exbins).astype(int))], axis=1)
	# H_extended[:,grouping_Egamma:-1:grouping_Egamma] /= 2
	# H_extended[:,grouping_Egamma+1:-2:grouping_Egamma] /= 2
	# H_extended = H[:,0:i_Egamma_max].repeat(N_Exbins).reshape(len(H[:,0]),N_Exbins,i_Egamma_max).sum(axis=2)/N_Exbins
	# plt.matshow(H_extended)
	# plt.show()
	# H_compressed = H[:,0:i_Egamma_max+ N_Exbins - i_Egamma_max%N_Exbins].reshape(N_Exbins, N_Exbins, grouping_Egamma).sum(axis=2)
	# plt.matshow(H_compressed)
	# plt.show()
	# H_compressed = rebin(H[:,0:i_Egamma_max], N_Exbins, 1)
	# plt.matshow(H_compressed)
	# plt.show()

	# H_compressed_extended = H_extended.reshape(N_Exbins, N_Exbins, grouping_Egamma).sum(axis=2)
	# plt.matshow(H_compressed_extended)
	# plt.show()

	# sys.exit(0)

	# Declare variables which will define the limits for the diff spectrum colorbar (for plotting purposes)
	vmin_spec = -200
	vmax_spec = 200
	vmin_diff = -100
	vmax_diff = 100

	# Perform the iterative subtraction:
	for iteration in range(N_iterations):
	# convergence_criterion = 1
	# max_diff = 100
	# while max_diff > convergence_criterion:
		# Store H from previous iteration to compare at the end
		H_old = H
		# Compress the H matrix along gamma axis to facilitate conversion to excitation energy
		# H_compressed = H[:,0:i_Egamma_max].reshape(N_Exbins, N_Exbins, grouping_Egamma).sum(axis=2)
		H_compressed = rebin(H[:,0:i_Egamma_max], N_Exbins, rebin_axis=1)

		# plt.pcolormesh(Egamma_range_compressed, Ex_range, H_compressed)
		# plt.show()

		# Convert first-generation spectra H into weights W
		W = np.zeros((N_Exbins, N_Exbins))
		for i in range(0,N_Exbins):
			# print H_compressed[i,i:0:-1].shape
			W[i,0:i] = H_compressed[i,i:0:-1]
		# plt.matshow(W, origin='lower', vmin=W.min(), vmax=W.max())
		# plt.colorbar()
		# plt.title('Before')
		# plt.show()
		# Remove negative weights
		W[W<0] = 0
		# Normalize each Ex channel to unity
		# W = np.where(np.invert(np.isnan(W/W.sum(axis=1).astype(float))),  W/W.sum(axis=1).astype(float), 0)
		# Remove Inf and NaN
		W = div0(W, W.sum(axis=1).reshape(N_Exbins,1))
		# W = np.nan_to_num(W) 
		# plt.matshow(W, origin='lower', vmin=W.min(), vmax=W.max())
		# plt.colorbar()
		# plt.title('After')
		# plt.show()

		# sys.exit(0)

		# print "W = "
		# print W
		# print "matrix_ex_compressed = "
		# print matrix_ex_compressed
		# print "product ="
		# plt.matshow(np.dot(W, matrix_ex_compressed), origin='lower', norm=LogNorm())
		# plt.show()

		# Calculate product of normalization matrix, weight matrix and raw count matrix
		G = np.dot( (normalization_matrix * W), matrix_ex_compressed) # Matrix of weighted sum of spectra below
		
		# Apply area correction
		if apply_area_correction:
			# Setup meshgrids for making boolean indexing arrays
			# Egamma_mesh_compressed, Ex_mesh_compressed = np.meshgrid(Egamma_range_compressed, Ex_range)
			# Egamma_max = Ex_range + dE_gamma # Maximal Egamma value for each Ex bin
			# Egamma_max_grid_compressed = np.meshgrid(np.ones(N_Exbins), Egamma_max)[1]
			# print "Egamma_mesh_compressed, Egamma_max, Egamma_max_grid"
			# print Egamma_mesh_compressed.shape, Egamma_max.shape, Egamma_max_grid.shape
			if statistical_or_total == 1:
				# Statistical multiplicity calculation (i.e. trying to use statistical/continuum region only)
				# slide_compressed = np.minimum( np.maximum(ThresRatio*Ex_mesh_compressed, ThresTot), ThresSta ) # The sliding lower limit for Egamma integral - sliding between ThresTot and ThresSta.
				# print "slide_compressed"
				# print slide_compressed.shape
				# plt.figure(5)
				# plt.plot(slide[:,0])
				# plt.show()
				# sys.exit(0)
				# good_indices_G = np.where(np.logical_and(slide_compressed < Egamma_mesh_compressed, Egamma_mesh_compressed < Egamma_max_grid) , True, False)
				G_area = np.where(np.logical_and(slide < Egamma_mesh, Egamma_mesh < Egamma_max_grid), G, 0).sum(axis=1)
			elif statistical_or_total == 2:
				# Total multiplicity calculation
				# good_indices_G = np.where(Egamma_mesh_compressed < Egamma_max_grid, True, False)
				G_area = np.where(Egamma_mesh_compressed < Egamma_max_grid, G, 0).sum(axis=1)
			# G_area = np.where(good_indices_G, G, 0).sum(axis=1)
			# print "print G_area.shape"
			# print G_area.shape
			# print "print G_area"
			# print G_area
			alpha = np.where(G_area > 0, (1 - div0(1,multiplicity)) * div0( area_matrix_ex_compressed_cut, G_area ), 1)
			alpha[alpha < 0.85] = 0.85
			alpha[alpha > 1.15] = 1.15
			# print "alpha.shape"
			# print alpha.shape
			# print "alpha"
			# print alpha
		else:
			alpha = np.ones(N_Exbins)


		# The actual subtraction
		H = matrix_ex_compressed - alpha.reshape((len(alpha), 1))*G
		# print H.shape
		# Plotting:
		# vmin_diff = (H-H_old).min()
		# vmax_diff = (H-H_old).max()
		# vmin_spec = H.min()
		# vmax_spec = H.max()

		# plt.figure(10)
		# plt.subplot(1,2,1)
		# plt.title('First gen spectrum, current')
		# plt.pcolormesh(Egamma_range, Ex_range, H, norm=LogNorm(vmin=0.01, vmax=vmax_spec))
		# plt.colorbar()
		# plt.subplot(1,2,2)
		

		# plt.title('Diff with previous')
		# plt.pcolormesh(Egamma_range, Ex_range, H-H_old, vmin=vmin_diff, vmax=vmax_diff)
		# plt.colorbar()
		# plt.show()



		# Check convergence
		max_diff = np.max(np.power(H-H_old,2))
		# print max_diff

	
	# Remove negative counts
	H[H<0] = 0
	# Return
	return H, H-H_old, Egamma_range, Ex_range




def rebin(array, N_final, rebin_axis=0):
	# Function to rebin an M-dimensional array either to larger or smaller binsize.
	# Rebinning is done with simple proportionality. E.g. for down-scaling rebinning (N_final < N_initial): 
	# if a bin in the original spacing ends up between two bins in the reduced spacing, 
	# then the counts of that bin is split proportionally between adjacent bins in the 
	# rebinned array. 
	# Upward binning (N_final > N_initial) is done in the same way, dividing the content of bins
	# equally among adjacent bins.

	# Technically it's done by repeating each element of array N_final times and dividing by N_final to 
	# preserve total number of counts, then reshaping the array from M dimensions to M+1 before flattening 
	# along the new dimension of length N_initial, resulting in an array of the desired dimensionality.
	indices = np.insert(array.shape, rebin_axis, N_final) # Indices to reshape 
	return array.repeat(N_final, axis=rebin_axis).reshape(indices).sum(axis=(rebin_axis+1))/float(N_final)





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

	if isinstance(array, tuple): # Check if input array is actually a tuple, which may happen if rebin_and_shift() is called several times nested for different axes.
		array = array[0]

	
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





def div0( a, b ):
    """ division function designed to ignore / 0, i.e. div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c








def rhosigchi(fgmat, fgvar, Egamma_range, Ex_range, N, method="BFGS"):
	##########################################################################
	# Test rhosigchi fitting:
	##########################################################################
	
	# filename_fg = "firstgen_test1.m"
	# filename_fg = "test.m"
	# fgmat, a, x_array, y_array = read_mama(filename_fg)
	# fgmat[fgmat < 0] = 0
	# fgmat = firstgen_matrix
	
	# Plot to check
	# plt.pcolormesh(x_array, y_array, fgmat, norm=LogNorm(vmin=0.001, vmax=fgmat.max()), cmap='gist_rainbow_r')
	# # plt.matshow(matrix)
	# # plt.pcolormesh(matrix, norm=LogNorm(vmin=0.1, vmax=matrix.max()), cmap='gist_rainbow_r')
	# plt.colorbar()
	# plt.xlabel('$E_\gamma$ [keV]', fontsize=14)
	# plt.ylabel('$E_x [keV]$', fontsize=14)
	# # plt.ylim([0,12])
	# # plt.xlim([0,12])
	# plt.show()
	
	# sys.exit(0)
	
	# print a
	# print x_array, y_array
	
	# Egamma_mesh, Ex_mesh = np.meshgrid(x_array, y_array)
	# TODO: Move these settings into arguments / read from input file
	Ex_high = 7500
	Ex_low = 4000
	Egamma_low = 1500
	Egamma_padding = 900 # We allow gamma energies somewhat higher than Ex_max because of detector uncertainty

	
	# Calculate maximal gamma energy we include, given by Ex_max + padding
	Egamma_max = Ex_range[-1] + Egamma_padding
	if Egamma_max > Egamma_range[-1]:
		Egamma_max = Egamma_range[-1]
	# Take out the slice along gamma axis corresponding to these energy limits (also for variance matrix)
	fgmat_padded, Egamma_range_sliced = slice_matrix_simple(fgmat, Egamma_range, [0,Egamma_max], axis=1)
	fgvar_padded, tmp = slice_matrix_simple(fgvar, Egamma_range, [0,Egamma_max], axis=1)
	plt.figure(0)
	plt.pcolormesh(Egamma_range_sliced, Ex_range, fgmat_padded)
	plt.show()
	# Add the same padding on top of Ex to get the same max energy
	N_Ex_padding_top = int(Egamma_padding/(Ex_range[1]-Ex_range[0])) # Number of additional Ex bins needed on top to reach same max energy
	# Also we need it below zero in both Ex and Egamma direction to facilitate transformation to Ef-Egamma coordinates
	N_Ex_padding_bottom = int((Egamma_padding+Ex_range[0])/(Ex_range[1]-Ex_range[0])) 								   # Number of additional Ex/Eg bins needed below zero,
	N_Eg_padding_bottom = int((Egamma_padding+Egamma_range_sliced[0])/(Egamma_range_sliced[1]-Egamma_range_sliced[0])) # taking into account that first Ex point may not be zero.
	# Apply the padding
	fgmat_padded = np.append(fgmat_padded, np.zeros((N_Ex_padding_top, fgmat_padded.shape[1])), axis=0) # Pad with zeros above in Ex
	fgmat_padded = np.append(np.zeros((N_Ex_padding_bottom, fgmat_padded.shape[1])), fgmat_padded, axis=0) # Pad with zeros below in Ex
	fgmat_padded = np.append(np.zeros((fgmat_padded.shape[0], N_Eg_padding_bottom)), fgmat_padded, axis=1) # Pad with zeros below in Eg
	# Same for fgvar
	fgvar_padded = np.append(fgvar_padded, np.zeros((N_Ex_padding_top, fgvar_padded.shape[1])), axis=0) # Pad with zeros above in Ex
	fgvar_padded = np.append(np.zeros((N_Ex_padding_bottom, fgvar_padded.shape[1])), fgvar_padded, axis=0) # Pad with zeros below in Ex
	fgvar_padded = np.append(np.zeros((fgvar_padded.shape[0], N_Eg_padding_bottom)), fgvar_padded, axis=1) # Pad with zeros below in Eg
	# Rebin to NxN
	fg_EiEg = rebin(rebin(fgmat_padded, N, rebin_axis=0), N, rebin_axis=1)
	fgv_EiEg = rebin(rebin(fgvar_padded, N, rebin_axis=0), N, rebin_axis=1)
	# Make the corresponding axis arrays of length N
	Ex_range_squared = np.linspace(Ex_range[0]-Egamma_padding, Ex_range[-1]+Egamma_padding, N)
	Eg_range_squared = np.linspace(Egamma_range_sliced[0]-Egamma_padding, Egamma_range_sliced[-1], N)

	# In case we want to apply some more restricting energy cuts than in the FG spectrum, they are applied here using a 2D meshgrid of energy ranges:
	Egammamesh, Exmesh = np.meshgrid(Eg_range_squared, Ex_range_squared)
	# Create a binary mask defining the area of the fg matrix that we are restricting to
	mask_EiEg = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), 1, 0) 
	# Also make the mask cut away gamma counts higher than Ex + E_gamma_padding
	mask_EiEg = np.where( Egammamesh < Exmesh + Egamma_padding, mask_EiEg, 0 ) 
	# Apply mask to firstgen matrix and variance matrix:
	fg_EiEG = mask_EiEg*fg_EiEg
	fgv_EiEg = mask_EiEg*fgv_EiEg
	# fgmat_cut = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), fgmat_padded, 0)
	# fgvar_cut = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), fgvar_padded, 0)
	plt.pcolormesh(fg_EiEg)
	
	plt.subplot(2,2,3)
	# Normalize the fg matrix and variance matrix (for each Ei bin, which is the way to do it to have normalized branching ratios)
	fg_EiEg = div0(fg_EiEg, fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1))
	fgv_EiEg = div0(fgv_EiEg, np.power(fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1),2)) # Normalize the variance accordingly, using that Var(aX) = a^2 Var(X).
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fg_EiEg)
	# # Need two different versions of matrices, one with Ei and one with Ef on y axis.
	fg_EfEg = EitoEf(fg_EiEg, Ex_range_squared)
	fgv_EfEg = EitoEf(fgv_EiEg, Ex_range_squared)
	# Also get the mask in EfEg coordinates
	mask_EfEg = EitoEf(mask_EiEg, Ex_range_squared)


	# Test applying smoothing to the variance matrix:
	plt.figure(5)
	plt.subplot(2,2,1)
	plt.title('Smoothing filter')
	print "fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	fgv_EiEg[fgv_EiEg < 4] = 4 # No errors less than sqrt(4)=2 (from rhosigchi.f)
	plt.pcolormesh(fgv_EiEg)
	from scipy.signal import savgol_filter
	fgv_EiEg = savgol_filter(fgv_EiEg, 5, 3, axis=1)
	print "smoothed fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.subplot(2,2,2)
	plt.pcolormesh(fgv_EiEg)
	plt.subplot(2,2,3)
	# Apply mask to smoothed fgv
	fgv_EiEg = mask_EiEg*fgv_EiEg
	print "smoothed and masked fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.pcolormesh(fgv_EiEg)
	# Make sure no negative variances (the fact that these even appear suggest that smoothing variances isn't very nice to do)
	fgv_EiEg[fgv_EiEg < 0] = 0
	print "smoothed, masked and removed negative fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.subplot(2,2,4)
	plt.pcolormesh(fgv_EiEg)
	plt.show()
	
	
	# fgmat_EiEgamma = np.flipud(fgmat_EiEgamma)
	# print "fgmat.sum(axis=1)"
	# print fgmat.sum(axis=1)
	# plt.pcolormesh(fgmat_EiEgamma)
	
	# plt.subplot(2,2,4)
	# # Normalize
	# fgmat_EiEgamma = div0(fgmat_EiEgamma, fgmat_EiEgamma.sum(axis=1).reshape(fgmat_EiEgamma.shape[0],1))
	# plt.pcolormesh(fgmat_EiEgamma)
	
	
	
	# plt.show()
	
	
	# sys.exit(0)
	
	
	
	def error(x):
		N = int(len(x)/2)
		T = x[0:N]
		rho = x[N:2*N]
		T2D, rho2D = np.meshgrid(T, rho, indexing='xy')
		fg_fit = T2D*rho2D
		chisquare = np.sum( div0(np.power(fg_fit-fg_EfEg, 2), fgv_EfEg) ) 
		# print chisquare
		return chisquare
	
	from scipy.optimize import minimize
	rho0 = np.ones(N)
	# T0 = np.linspace(0,1,N)
	# T0 = rebin(fg_Egamma.mean(axis=0), N, rebin_axis=0)
	T0 = fg_EfEg.sum(axis=0)
	x = np.append(T0, rho0)
	
	T = x[0:N]
	rho = x[N:2*N]
	T2D, rho2D = np.meshgrid(T, rho, indexing='xy')
	# print "T2D"
	# plt.matshow(T2D)
	# print "rho2D"
	# print rho2D
	# plt.plot(T0)
	# plt.show()
	# print x
	# print error(x)
	
	result = minimize(error, x, method=method, options={"maxiter":1e6,"maxfev":1e6})
	print "Rhosigchi minimization finished. Obtained convergence (True/False):", result.success
	if result.success:
		Nfree = fgmat.shape[0]*fgmat.shape[1]/2 - len(rho0) - len(T0)
		print "Minimal Chi-square:", 1/Nfree * result.fun, ", number of iterations:", result.nit, "degrees of freedom:", Nfree
	elif not result.success:
		print result.message
	
	if result.success:
		T_fit = result.x[0:N]
		rho_fit = result.x[N:2*N]
		# print T_fit, rho_fit
		
		plt.figure(20)
		plt.subplot(2,1,1)
		Egamma_plot_range = np.linspace(1500, 8400, N)
		Ex_plot_range = np.linspace(0, 7500, N)
		alpha = 0.0018
		# alpha = 0
		# plt.plot(Egamma_plot_range, T_fit/Egamma_plot_range**3 * np.exp(alpha*Egamma_plot_range), 'o')
		plt.plot(Egamma_plot_range, T_fit, 'o')
		plt.yscale('log')
		plt.subplot(2,1,2)
		# plt.plot(Ex_plot_range, rho_fit* np.exp(alpha*Ex_plot_range), 'o')
		plt.plot(Ex_plot_range, rho_fit, 'o')
		plt.yscale('log')
		plt.show()
	




























def rhosigchi2(fgmat, fgvar, Egamma_range, Ex_range, dE_gamma, N):
	# Test rhosigchi fitting. This version uses Schiller's (NIM, 2000) iterative method.

	# TODO: Move these settings into arguments / read from input file
	Ex_high = 7500
	Ex_low = 4000
	Egamma_low = 1500
	Egamma_padding = 900 # We allow gamma energies somewhat higher than Ex_max because of detector uncertainty
	
	# plt.figure(0)
	# plt.subplot(2,2,1)
	print "Calibration coefficients:"
	print "bx =", Ex_range[0], "ax =", Ex_range[1]-Ex_range[0]
	print "bg =", Egamma_range[0], "ag =", Egamma_range[1]-Egamma_range[0]

	# Calculate maximal gamma energy we include, given by Ex_max + padding
	Egamma_max = Ex_range[-1] + Egamma_padding
	if Egamma_max > Egamma_range[-1]:
		Egamma_max = Egamma_range[-1]
	# Take out the slice along gamma axis corresponding to these energy limits (also for variance matrix)
	fgmat_padded, Egamma_range_sliced = slice_matrix_simple(fgmat, Egamma_range, [0,Egamma_max], axis=1)
	fgvar_padded, tmp = slice_matrix_simple(fgvar, Egamma_range, [0,Egamma_max], axis=1)
	plt.figure(0)
	plt.pcolormesh(Egamma_range_sliced, Ex_range, fgmat_padded)
	plt.show()
	# Add the same padding on top of Ex to get the same max energy
	N_Ex_padding_top = int(Egamma_padding/(Ex_range[1]-Ex_range[0])) # Number of additional Ex bins needed on top to reach same max energy
	# Also we need it below zero in both Ex and Egamma direction to facilitate transformation to Ef-Egamma coordinates
	N_Ex_padding_bottom = int((Egamma_padding+Ex_range[0])/(Ex_range[1]-Ex_range[0])) 								   # Number of additional Ex/Eg bins needed below zero,
	N_Eg_padding_bottom = int((Egamma_padding+Egamma_range_sliced[0])/(Egamma_range_sliced[1]-Egamma_range_sliced[0])) # taking into account that first Ex point may not be zero.
	# Apply the padding
	fgmat_padded = np.append(fgmat_padded, np.zeros((N_Ex_padding_top, fgmat_padded.shape[1])), axis=0) # Pad with zeros above in Ex
	fgmat_padded = np.append(np.zeros((N_Ex_padding_bottom, fgmat_padded.shape[1])), fgmat_padded, axis=0) # Pad with zeros below in Ex
	fgmat_padded = np.append(np.zeros((fgmat_padded.shape[0], N_Eg_padding_bottom)), fgmat_padded, axis=1) # Pad with zeros below in Eg
	# Same for fgvar
	fgvar_padded = np.append(fgvar_padded, np.zeros((N_Ex_padding_top, fgvar_padded.shape[1])), axis=0) # Pad with zeros above in Ex
	fgvar_padded = np.append(np.zeros((N_Ex_padding_bottom, fgvar_padded.shape[1])), fgvar_padded, axis=0) # Pad with zeros below in Ex
	fgvar_padded = np.append(np.zeros((fgvar_padded.shape[0], N_Eg_padding_bottom)), fgvar_padded, axis=1) # Pad with zeros below in Eg
	# Rebin to NxN
	fg_EiEg = rebin(rebin(fgmat_padded, N, rebin_axis=0), N, rebin_axis=1)
	fgv_EiEg = rebin(rebin(fgvar_padded, N, rebin_axis=0), N, rebin_axis=1)
	# Make the corresponding axis arrays of length N
	Ex_range_squared = np.linspace(Ex_range[0]-Egamma_padding, Ex_range[-1]+Egamma_padding, N)
	Eg_range_squared = np.linspace(Egamma_range_sliced[0]-Egamma_padding, Egamma_range_sliced[-1], N)

	# In case we want to apply some more restricting energy cuts than in the FG spectrum, they are applied here using a 2D meshgrid of energy ranges:
	Egammamesh, Exmesh = np.meshgrid(Eg_range_squared, Ex_range_squared)
	# Create a binary mask defining the area of the fg matrix that we are restricting to
	mask_EiEg = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), 1, 0) 
	# Also make the mask cut away gamma counts higher than Ex + E_gamma_padding
	mask_EiEg = np.where( Egammamesh < Exmesh + Egamma_padding, mask_EiEg, 0 ) 
	# Apply mask to firstgen matrix and variance matrix:
	fg_EiEg = mask_EiEg*fg_EiEg
	fgv_EiEg = mask_EiEg*fgv_EiEg
	# fgmat_cut = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), fgmat_padded, 0)
	# fgvar_cut = np.where(np.logical_and(Exmesh > Ex_low, np.logical_and(Exmesh < Ex_high, Egammamesh > Egamma_low)), fgvar_padded, 0)
	plt.pcolormesh(fg_EiEg)
	
	plt.subplot(2,2,3)
	# Normalize the fg matrix and variance matrix (for each Ei bin, which is the way to do it to have normalized branching ratios)
	fg_EiEg = div0(fg_EiEg, fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1))
	fgv_EiEg = div0(fgv_EiEg, np.power(fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1),2)) # Normalize the variance accordingly, using that Var(aX) = a^2 Var(X).
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fg_EiEg)
	# # Need two different versions of matrices, one with Ei and one with Ef on y axis.
	fg_EfEg = EitoEf(fg_EiEg, Ex_range_squared)
	fgv_EfEg = EitoEf(fgv_EiEg, Ex_range_squared)
	# Also get the mask in EfEg coordinates
	mask_EfEg = EitoEf(mask_EiEg, Ex_range_squared)
	

	# fgmat_EiEgamma = rebin(fg_EiEg, fgvar_padded.shape[0], rebin_axis=1) # TODO: Think carefully about the math in rebinning a variance matrix. Does it make sense? Probably OK as long as it's the variance and not stdev.
	# fgvar_EiEgamma = fgv_EiEg
	# fgmat_EfEgamma = np.zeros(fg_EiEg.shape)
	# fgvar_EfEgamma = np.zeros(fgv_EiEg.shape)
	# for i in range(fg_EiEg.shape[0]):
	# 	fgmat_EfEgamma[i,0:fg_EiEg.shape[0]-i] = fg_EiEg.diagonal(-i)
	# 	fgvar_EfEgamma[i,0:fg_EiEg.shape[0]-i] = fgv_EiEg.diagonal(-i)
	# fgmat_EfEgamma = EitoEf(fgmat_EiEgamma, Ex_range_padded)
	# fgvar_EfEgamma = EitoEf(fgmat_EiEgamma, Ex_range_padded)
	plt.subplot(2,2,4)
	plt.pcolormesh(fg_EfEg)
	# plt.show()

	plt.figure(6)
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, mask_EiEg)
	plt.title('mask_EiEg')
	plt.show()

	# Map back to EiEg to check if my inverse map is right:
	# control_EiEgamma = np.zeros(fg_EiEg.shape)
	# for i in range(fg_EiEg.shape[0]):
	# 	control_EiEgamma[fg_EiEg.shape[0]-i-1,0:fg_EiEg.shape[0]-i] = np.flipud(fgmat_EfEgamma).diagonal(-i)
	# plt.figure(11)
	# plt.pcolormesh(control_EiEgamma)
	# plt.show()

	# sys.exit(0)

	# Chop away zeros. Is this really necessary? Test without. Should in any case be done generally, using the energy limits set earlier.
	# fgmat_EfEgamma_cut = fgmat_EiEgamma[0:150,55:205]
	# fgvar_EfEgamma_cut = fgvar_EiEgamma[0:150,55:205]

	# plt.matshow(fgmat_EiEgamma_cut)
	# fg_EiEg = rebin(rebin(fgmat_EiEgamma, N, rebin_axis=0), N, rebin_axis = 1)
	# fgv_EiEg = rebin(rebin(fgvar_EiEgamma, N, rebin_axis=0), N, rebin_axis = 1)
	# Ex_range_padded_rebinned = np.linspace(Ex_range_padded[0], Ex_range_padded[-1], N)
	# Eg_range_padded_rebinned = np.linspace(Eg_range_padded[0], Eg_range_padded[-1], N)
	# mask_EiEg = np.where(rebin(rebin(mask_EiEg, N, rebin_axis=0), N, rebin_axis = 1) > 0, 1, 0) # Rebin the masking array, resetting to binary values in an inclusive way (i.e. edge values are upscaled to 1 rather than down to 0)
	# mask_EfEg = EitoEf(mask_EiEg)
	# plt.pcolormesh(Eg_range_padded_rebinned, Ex_range_padded_rebinned, mask_EfEg)
	plt.show()


	# sys.exit(0)

	# Test applying smoothing to the variance matrix:
	plt.figure(5)
	plt.subplot(2,2,1)
	plt.title('Smoothing filter')
	print "fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	fgv_EiEg[fgv_EiEg < 4] = 4 # No errors less than sqrt(4)=2 (from rhosigchi.f)
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fgv_EiEg)
	from scipy.signal import savgol_filter
	fgv_EiEg = savgol_filter(fgv_EiEg, 21, 3, axis=1)
	print "smoothed fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.subplot(2,2,2)
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fgv_EiEg)
	plt.subplot(2,2,3)
	# Apply mask to smoothed fgv
	fgv_EiEg = mask_EiEg*fgv_EiEg
	print "smoothed and masked fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fgv_EiEg)
	# Make sure no negative variances (the fact that these even appear suggest that smoothing variances isn't very nice to do)
	fgv_EiEg[fgv_EiEg < 0] = 0
	print "smoothed, masked and removed negative fgv: mean =", fgv_EiEg.mean(), "min =", fgv_EiEg.min(), "max =", fgv_EiEg.max()
	plt.subplot(2,2,4)
	plt.pcolormesh(Eg_range_squared, Ex_range_squared, fgv_EiEg)
	plt.show()


	# # HACK TEST: Set fgv = fg. Except for scaling by multiplicities, this is what Schiller does. Seems to give crazy results in this case, must be something wrong!
	# fgv_EiEg = 2*np.random.normal(loc=1, scale=0.01, size=fg_EiEg.shape)*fg_EiEg


	# Get EfEg arranged versions of matrices
	# fg_EfEg =  EitoEf(fg_EiEg, Ex_range_padded)
	# fgv_EfEg = EitoEf(fgv_EiEg, Ex_range_padded)


	# plt.figure(1)
	# plt.subplot(2,1,1)
	# plt.pcolormesh(fg_EiEg)
	# plt.subplot(2,1,2)
	# plt.pcolormesh(fg_EfEg)
	# plt.show()

	# plt.pcolormesh(fgv_EfEg)
	# plt.show()

	# sys.exit(0)


	# Make initial F and rho
	rho = np.ones(N)
	# F = rebin(fg_EfEg.mean(axis=0), N, rebin_axis=0)
	F = fg_EiEg.mean(axis=0) # CHECK THAT THIS IS PROPERLY NORMALIZED
	# plt.figure(2)
	# plt.plot(F)
	# plt.show()

	# Now start iterating.
	N_iterations = 2
	for iteration in range(N_iterations):

		# Plot current estimate of F and rho
		plt.figure(20)
		plt.title('F (top) and rho (bottom). iteration %d'%iteration)
		plt.subplot(2,1,1)
		plt.plot(F)
		plt.yscale('log')
		plt.subplot(2,1,2)
		plt.plot(rho)
		plt.yscale('log')
		plt.show()



		# Make meshgrid of F and rho
		rho2D_EfEg, F2D_EfEg = np.meshgrid(rho, F, indexing='ij')
		# rho2D_EfEg = rho2D_EfEg.T # Transpose to make it right?

		# Also need them as functions of Ei:
		# F2D_EiEg = np.zeros(F2D_EfEg.shape)
		# rho2D_EiEg = np.zeros(rho2D_EfEg.shape)
		# for i in range(F2D_EfEg.shape[0]):
		# 	F2D_EiEg[F2D_EfEg.shape[0]-i-1,0:F2D_EfEg.shape[0]-i] = np.flipud(F2D_EfEg).diagonal(-i)
		# 	rho2D_EiEg[F2D_EfEg.shape[0]-i-1,0:F2D_EfEg.shape[0]-i] = np.flipud(rho2D_EfEg).diagonal(-i)
		# 	print "check"
		F2D_EiEg = EftoEi(F2D_EfEg, Ex_range_squared)
		rho2D_EiEg = EftoEi(rho2D_EfEg, Ex_range_squared)

		# # Plot to inspect the F in Ef and Ei coordinates:
		# plt.figure(3)
		# plt.subplot(2,1,1)
		# plt.pcolormesh(F2D_EfEg)
		# plt.subplot(2,1,2)
		# plt.pcolormesh(F2D_EiEg)
		# plt.show()	

		# Plot the 2D product of F and rho in EiEg coordinates (the "theoretical first-generation spectrum"):
		fgtheo_EiEg = F2D_EiEg*rho2D_EiEg
		plt.figure(4)
		plt.title('Theoretical first generation spectrum (product of F and rho)')
		plt.subplot(2,2,1)
		plt.title('F2D_EiEg')
		plt.pcolormesh(Eg_range_squared, Ex_range_squared, F2D_EiEg)
		plt.subplot(2,2,2)
		plt.title('rho2D_EiEg')
		plt.pcolormesh(Eg_range_squared, Ex_range_squared, rho2D_EiEg)
		plt.subplot(2,2,3)
		plt.title('fgtheo_EiEg')
		plt.pcolormesh(Eg_range_squared, Ex_range_squared, fgtheo_EiEg*mask_EiEg)
		plt.subplot(2,2,4)
		plt.title('fg_EiEg')
		plt.pcolormesh(Eg_range_squared, Ex_range_squared, fg_EiEg*mask_EiEg)
		plt.show()
	
		# Now calculate all the helping functions(vectors/matrices) according to Schiller (NIM, 2000), with the same naming:
		# s = ( F2D_EiEg * rho2D_EfEg *mask_EiEg).sum(axis=1) # s(Ei) # CHANGED 20161020 rho was EiEg should be EfEg? Also tried changing F to EfEg but unsure. 20161022 changed back for further testing, also multiplied with mask. Now also trying to use different for the two (Ei for F and Ef for rho)
		# plt.figure()
		# plt.subplot(1,2,1)
		# plt.title('F2D_EiEg*rho2D_EfEg ')
		# plt.pcolormesh(Eg_range_squared, Ex_range_squared, F2D_EiEg * rho2D_EfEg * mask_EiEg)
		# plt.subplot(1,2,2)
		# plt.title('F2D*rho2D EiEg both')
		# plt.pcolormesh(Eg_range_squared, Ex_range_squared, F2D_EiEg * rho2D_EiEg * mask_EiEg)

		# 20161023: Trying to manually calculate the s, a and b with for loops:


		plt.figure()
		plt.plot(s, label='EiEg')
		plt.plot(( F2D_EfEg * rho2D_EfEg * mask_EfEg).sum(axis=1), label='EfEg')
		plt.title('s comparison, iteration %d'%iteration)
		plt.legend()
		plt.show()
		# plt.matshow(np.log(np.power( div0( F2D_EiEg*rho2D_EiEg, fgv_EiEg ), 2)))
		# plt.show()
		a = ( div0( np.power(F2D_EiEg*rho2D_EfEg, 2), fgv_EiEg )  *mask_EiEg).sum(axis=1) # a(Ei) # CHANGED 20161020 rho was EiEg should be EfEg? Also tried changing F to EfEg but unsure. 20161022 changed back for further testing, also multiplied with mask.
		b = ( div0( F2D_EiEg*rho2D_EfEg*fg_EiEg, fgv_EiEg )  *mask_EiEg).sum(axis=1) # b(Ei) # CHANGED 20161020 rho was EiEg should be EfEg? Also tried changing F to EfEg but unsure. 20161022 changed back for further testing, also multiplied with mask.
		# plt.matshow(np.log(div0( F2D_EiEg*rho2D_EiEg*fg_EiEg, np.power(fgv_EiEg, 2) )))
		# plt.show()
		# Need matrix versions for multiplication. Each column is a repeat of the s/a/b vector, respectively
		s2D = s.repeat(N).reshape((N,N))
		# plt.matshow(s2D)
		# plt.show()
		a2D = a.repeat(N).reshape((N,N))
		# plt.matshow(a2D)
		# plt.show()
		b2D = b.repeat(N).reshape((N,N))
		# plt.matshow(np.log(b2D))
		# plt.show()
		psi = div0( np.ones(fg_EiEg.shape), np.power( s2D, 2 )*fgv_EiEg )  # psi(Ei,Eg)
		psi = mask_EiEg * psi
		# plt.matshow(np.log(psi))
		# plt.title('log(psi)')
		# plt.show()
		print "Checking numerical stability: a.mean() =", a.mean(), "b.mean() =", b.mean(), "s.mean() =", s.mean()
		phi = div0( a2D, np.power(s2D, 3) ) - div0( b2D, np.power(s2D, 2) ) + div0( fg_EiEg, s2D*fgv_EiEg ) # phi(Ei,Eg)
		print "a/s^3.mean() =", div0( a2D, np.power(s2D, 3) ).mean(), "b/s^2.mean() =", div0( b2D, np.power(s2D, 2) ).mean(), "fg/s*fgv^2.mean() =", div0( fg_EiEg, s2D*fgv_EiEg ).mean()
		phi = mask_EiEg * phi # Apply mask
		# plt.matshow(np.log(phi))
		# plt.title('log(phi)')
		# plt.show()

		# Store old:
		F_old = F       
		rho_old = rho   

		# Calculate updated F and rho:

		# HACK: Multiplying by 1e2 to counter the unwanted behaviour where the function keeps decreasing through iterations. But it probably points to some more serious error somewhere.
		F =  div0( (rho2D_EfEg*phi).sum(axis=0), (np.power(rho2D_EfEg,2)*psi).sum(axis=0) ) # TODO: Could this be formulated as a matrix*vector operation?
		# Need to do some magic for rho: The relevant slices of phi and psi are diagonal Ei = Egamma + Ef. Mind-twisting! But I think it's simply the EfEg arranged versions of phi and psi.
		phiEfEg = EitoEf(phi, Ex_range_squared)
		psiEfEg = EitoEf(psi, Ex_range_squared)
		rho =  div0( np.dot(phiEfEg, F_old), np.dot(psiEfEg, np.power(F_old, 2)) )
		print "F.shape =", F.shape
		print "rho.shape =", rho.shape

		# Apply a test to ensure that rho and F don't change by more than 30 percent in either direction since last iteration
		# P = 0.7
		# if iteration > 0: # Don't do it for the first one
		# 	F = np.where(F_old/F > (1+P), F_old*(1+P), F)
		# 	F = np.where(F/F_old < 1/(1+P), F_old/(1+P), F)
		# 	rho = np.where(rho_old/rho > (1+P), rho_old*(1+P), rho)
		# 	rho = np.where(rho/rho_old < 1/(1+P), rho_old/(1+P), rho)






# def EitoEf(matrix_EiEg, Ex_range):
# 	matrix_EfEg = np.zeros(matrix_EiEg.shape)
# 	for i in range(matrix_EiEg.shape[0]):
# 		matrix_EfEg[i,0:matrix_EiEg.shape[0]-i] = matrix_EiEg.diagonal(-i)
# 	return matrix_EfEg

def EitoEf(matrix_EiEg, Ex_range):
	# Find out which index along Ex is the zero energy
	i_Ex_zero = np.where(Ex_range > 0)[0][0]
	if np.abs(Ex_range[i_Ex_zero]) > np.abs(Ex_range[i_Ex_zero-1]):
		i_Ex_zero -= 1

	# print "i_Ex_zero =", i_Ex_zero, "Ex_range[i_Ex_zero] =", Ex_range[i_Ex_zero]

	matrix_EfEg = np.zeros(matrix_EiEg.shape)
	# All the sub-diagonals from the main diagonal (Ei=Eg) and upwards in Ei > Eg direction are filled into the rows starting from i_Ex_zero:
	for i in range(0, matrix_EiEg.shape[0]-i_Ex_zero):
		matrix_EfEg[i+i_Ex_zero,0:len(matrix_EiEg.diagonal(-i))] = matrix_EiEg.diagonal(-i)
		# print matrix_EiEg.diagonal(-i)
	for i in range(1, i_Ex_zero+1):
		matrix_EfEg[i_Ex_zero-i,i:len(matrix_EiEg.diagonal(i))+i] = matrix_EiEg.diagonal(i)
	# print matrix_EiEg.shape
	# for i in range(0, matrix_EiEg.shape[0]): # i is the coordinate of the EfEg matrix row
	# 	matrix_EfEg[i,0:len(matrix_EiEg.diagonal(-i+i_Ex_zero))] = matrix_EiEg.diagonal(-i+i_Ex_zero)
		# print len(matrix_EiEg.diagonal(-i + i_Ex_zero))
	

	return matrix_EfEg

def EftoEi(matrix_EfEg, Ex_range):
	 # TODO test if this is proper inversion of EitoEf
	i_Ex_zero = np.where(Ex_range > 0)[0][0]
	if np.abs(Ex_range[i_Ex_zero]) > np.abs(Ex_range[i_Ex_zero-1]):
		i_Ex_zero -= 1
	# print "i_Ex_zero =", i_Ex_zero, "Ex_range[i_Ex_zero] =", Ex_range[i_Ex_zero]

	matrix_EiEg = np.zeros(matrix_EfEg.shape)
	for i in range(0, matrix_EfEg.shape[0]-i_Ex_zero):
		matrix_EiEg[matrix_EfEg.shape[0]-i-i_Ex_zero-1,0:len(np.flipud(matrix_EfEg).diagonal(-i))] = np.flipud(matrix_EfEg).diagonal(-i)
	for i in range(0, i_Ex_zero):
		matrix_EiEg[matrix_EfEg.shape[0]-i-1,matrix_EfEg.shape[0]-len(np.flipud(matrix_EfEg).diagonal(i_Ex_zero-i)):matrix_EfEg.shape[0]] = np.flipud(matrix_EfEg).diagonal(i_Ex_zero-i)
		# print i, -matrix_EfEg.shape[0]+i_Ex_zero+i, np.flipud(matrix_EfEg).diagonal(-matrix_EfEg.shape[0]+i+i_Ex_zero)
	# TODO: THe below probably needs to be included to account for Eg > Ei events. But must figure out how to do it. Use a test plot.
	# for i in range(1, i_Ex_zero):
	# 	matrix_EiEg[i,0:len(np.flipud(matrix_EfEg).diagonal(i))] = np.flipud(matrix_EfEg).diagonal(i)

	return matrix_EiEg

def slice_matrix_simple(matrix, E_range, E_limits, axis):
	# Returns a subslice of the matrix along specified axis, but does no rebinning. That means the slicing is only as precise as the bin spacing, 
	# since the nearest bin is used for the cut.
	Elow = E_limits[0]
	Ehigh = E_limits[1]
	index_Elow = np.where(E_range > Elow)[0][0]
	# print "E_range[index_Elow] first try =",E_range[index_Elow]
	if np.abs(Elow - E_range[index_Elow]) > np.abs(Elow - E_range[index_Elow-1]):
		index_Elow -= 1
	# print "E_range[index_Elow] second try =",E_range[index_Elow]
	index_Ehigh = np.where(E_range > Ehigh)[0][0]
	# print "E_range[index_Ehigh] first try =",E_range[index_Ehigh]
	if np.abs(Ehigh - E_range[index_Ehigh]) > np.abs(Ehigh - E_range[index_Ehigh-1]):
		index_Ehigh -= 1
	# print "E_range[index_Ehigh] second try =",E_range[index_Ehigh]
	# Return subslice of matrix along with corresponding subslice of axis range array E_range
	return np.split(matrix, [index_Elow, index_Ehigh+1], axis=axis)[1], E_range[index_Elow:index_Ehigh+1] # The +1 index is needed to actually include the final element that we want.

def slice_matrix(matrix, E_ranges, E_limits, N_finals):
	###########################################################
	# NB! This works, but it EATS memory. Use with caution!
	###########################################################
	# matrix is an N-dimensional histogram. 
	# E_ranges is a list of the arrays giving the center-bin energies corresponding to the axes: E_ranges = [E_range0, E_range1, ...] for axis 0, 1, ...
	# E_limits is a list of axis limits for the slice to be taken out of the matrix, on the form
	# E_limits = [E_low0, E_high0, E_low1, E_high1, ...], where 0, 1 etc. denote the axis along which the limits apply
	# Similarily N_finals is a list of the desired number of bins along corresponding axis after slice: N_finals = [N_final0, N_final1, ...]
	E_range_finals = []
	matrix_sliced = matrix
	for i in range(len(N_finals)):
		N_final = N_finals[i]
		Elow = E_limits[2*i]
		Ehigh = E_limits[2*i+1]
		E_range_original = E_ranges[i]
		# Make a more finely-spaced array, set up to have maximum 1 keV bin spacing.
		E_range_intermediate = np.linspace(E_range_original[0], E_range_original[-1], int(len(E_range_original)*(E_range_original[1]-E_range_original[0])*1))
		# Find the indices of E_range_intermediate that are closest to Elow and Ehigh, respectively
		index_Elow = np.where(E_range_intermediate > Elow)[0][0]
		if np.abs(Elow - E_range_intermediate[index_Elow]) > np.abs(Elow - E_range_intermediate[index_Elow-1]):
			index_Elow -= 1
		index_Ehigh = np.where(E_range_intermediate > Ehigh)[0][0]
		if np.abs(Ehigh - E_range_intermediate[index_Ehigh]) > np.abs(Ehigh - E_range_intermediate[index_Ehigh-1]):
			index_Ehigh -= 1
		# Make the axis energy bin array of the desired end product
		E_range_final = np.linspace(Elow, Ehigh, N_final)
		# Append the axis array to list of final axis arrays
		E_range_finals.append(E_range_final)
		# Rebin the matrix along current axis to get the same fine-grained spacing as for the intermediate array
		matrix_sliced = rebin(matrix_sliced, int(len(E_range_original)*(E_range_original[1]-E_range_original[0])*1), rebin_axis=i)
		# Take out the slice along current axis
		matrix_sliced = np.split(matrix_sliced, [index_Elow, index_Ehigh], axis=i)[1]
		# Rebin the slice to the desired spacing and name it matrix_sliced so it's picked up for the next axis round
		matrix_sliced = rebin(matrix_sliced, N_final, rebin_axis=i)
	return matrix_sliced, E_range_finals













# ===== OLD BELOW ======


# def first_generation_spectrum(matrix, limits, N_iterations = 20):
# 	Nx = len(matrix[0,:])
# 	Ny = len(matrix[:,0])

# 	Ngates= 0 # Number of excitation bins
# 	Iloop= 2	#  Iterate once more
# 	RDIM = 512			# Just to give flag that R(i,j) is used
# 	Iter       = 1			# Number of iterations
# 	ExH        = 8500.		# Highest excitation energy
# 	dE_gamma   = 500.         # Eg goes dE_gamma higher than Ex and lower than Ex=0 MeV, due to gamma-resolution
# 	ExpWeight  = 0			# Default not to weight with exp. spectra
# 	a          = 16.0			# Fermi level density parameter a (1/MeV)
# 	Nexp       = 4.2			# Exponent for gamma energy
# 	Norm       = 2			# Default normalization with singles method
# 	StaTot     = 1			# Choose statistical multiplicity
# 	ThresSta   = 430.			# with threshold > 430 keV
# 	ThresTot   = 200.			# Lower exp. threshold for gammas
# 	ThresRatio = 0.3			# Upper limit = AMIN1(Eg* ThresRatio,ThresSta)
# 	ExEntry0s  = 300.			# Average entry point in ground band for M_stat
# 	ExEntry0t  = 0.			# Average entry point in ground band for M_tot
# 	AreaCorr   = 0			# Correct for areas versus multiplicity
# 	ReadStatus = 0			# Has not read figegain.dat
# 	iAttnFac_active = 0		# By default, no active attenuation of direct decay

# 	# For testing purposes: Hard-coding the parameters from figegain.dat
# 	ExH = 7520
# 	Ax0 = 15.0
# 	Ax1 = 40.0
# 	Ay0 = 0.0
# 	Ay1 = 40.0
# 	Ngates = 189
# 	ExpWeight = 0
# 	A = 6.0
# 	Nexp = 4.19999981
# 	AxW0 = 0.0
# 	AxW1 = 0.0
# 	AyW0 = 0.0
# 	AyW1 = 0.0
# 	Norm = 2
# 	statistical_or_total = 1
# 	ThresSta = 300.0
# 	AreaCorr = 1
# 	Sing = np.array([0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000, 1000.00000])
# 	Mult = np.array([3.28284287, 3.33303022, 3.36863661, 3.38089156, 3.13876343, 3.20864224, 3.27202582, 3.18323016, 3.24140716, 3.18785071, 3.31042576, 3.18568373, 3.11389852, 3.17003512, 3.09723783, 3.15070581, 3.22956467, 3.20834374, 3.03619146, 2.99797702, 3.07418394, 2.95585322, 2.98966885, 3.02102375, 3.03368664, 2.95281339, 2.95928788, 2.85788822, 2.96565771, 2.89722013, 2.93754339, 2.86754513, 2.87389922, 2.94471216, 2.86285329, 2.88891768, 2.83815956, 2.89718437, 2.92867851, 2.88163972, 2.90771937, 2.88899922, 2.87331295, 2.80630994, 2.79502726, 2.85805321, 2.85137033, 2.96695566, 2.90483856, 2.80189657, 2.91725636, 2.85080647, 2.81791878, 2.73848462, 2.77616620, 2.75371099, 2.75522852, 2.64897299, 2.65639353, 2.75609660, 2.67403793, 2.63901448, 2.65125775, 2.63411927, 2.69224977, 2.66567993, 2.57934189, 2.62521434, 2.62242508, 2.61451650, 2.62248874, 2.50728250, 2.45229745, 2.36173773, 2.40619183, 2.40411615, 2.30891156, 2.31666422, 2.28118777, 2.38253140, 2.30362034, 2.33219433, 2.22109985, 2.27223325, 2.29964232, 2.22713542, 2.27693295, 2.30436993, 2.28863978, 2.31533909, 2.25484872, 2.29598260, 2.30992746, 2.25849462, 2.19287252, 2.16879559, 2.24799657, 2.18099260, 2.19578409, 2.20677543, 2.12944198, 2.17651463, 2.13324451, 2.09054422, 2.07191396, 2.02187300, 2.03894854, 2.06274152, 2.06018591, 2.04624677, 2.01378059, 1.93092191, 1.97711456, 1.94552124, 2.07609653, 2.01849103, 1.88140798, 1.94084358, 1.79872561, 1.91477597, 1.96998918, 1.88759911, 1.94787514, 1.90252268, 1.83544970, 1.78045821, 1.81808877, 1.77311456, 1.72606766, 1.74057031, 1.70578408, 1.65551829, 1.68457615, 1.69171178, 1.65942121, 1.59441555, 1.53900480, 1.64131331, 1.64370775, 1.60543394, 1.57258582, 1.62325704, 1.58347929, 1.55381584, 1.52986348, 1.47450376, 1.42342901, 1.49812043, 1.33530295, 1.36676729, 1.30415654, 1.38348877, 1.34240890, 1.32294381, 1.33545387, 1.37406445, 1.40528977, 1.37357938, 1.41531861, 1.38877332, 1.44060254, 1.30328941, 1.46965837, 1.40472591, 1.38291585, 1.35271358, 1.15609646, 1.14211118, 1.02502418, 1.00129557,0.963818729,0.909706414,0.959041893,0.950005472,0.981827021,0.937574327,0.897615790,0.839832962,0.755988777,0.728221178,0.711606681,0.661457121, 0.00000000,0.597014904,0.542372882,0.406779617, 0.00000000, 0.00000000, 0.00000000])
# 	AttnFac = np.array([0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000])
# 	ThresTot = 300.000000
# 	ThresRatio = 0.3000000
# 	ExH = 7520.00000
# 	ExEntry0s = 300.000000
# 	ExEntry0t = 0.00000000

# 	# Setup some parameters
# 	if StaTot == 1:
# 		ExEntry0 = ExEntry0s
# 	elif StaTot == 2:
# 		ExEntry0 = ExEntry0t
# 	else:
# 		print "Invalid value of StaTot"
# 		sys.exit(1)

# 	LS1 = int(((ThresSta-Ax0)/Ax1)+0.5)	 # Lower threshold for stat. gammas
# 	LT1 = int(((ThresTot-Ax0)/Ax1)+0.5)	 # Lower threshold for total gammas

# 	# Done with input parameters. Now allocate necessary matrices:
# 	# input_matrix, the total Ex-Eg matrix, is input as an argument to the function.
# 	R = np.zeros((2048,2048)) # Is this the response matrix? It's made by the Fermi Gas thing...


# 	# UPDATE: THIS IS REMOVED AND A SQUARE IS USED INSTEAD FOR THE FIRST ITERATION WEIGHT FUNCTION
# 	# # Fill response matrix according to Fermi gas model:
# 	IyH=max(int(((ExH-Ay0)/Ay1)+0.5), 0)
# 	IyL=max(int(((0-Ay0)/Ay1)+0.5), 0)
# 	IyH = min(IyH, 2048)
# 	IyL = min(IyL, 2048)
# 	IyStep=+1									
# 	if IyH>IyL:
# 		IyStep=-1
# 	Ngates=abs((IyH-IyL)+1)
# 	# Ex = Ay0 + Ay1*np.arange(IyH, IyL, IyStep) # Needs to be 2-D arrays (meshgrid?)
# 	# Egam = np.abs(Ay1)*np.arange(0,2048,1)	   # -"-
# 	# Egam[0] = np.abs(Ay1)*0.25 # To get something in ch 0
# 	# # print Ex, Egam
# 	# Ex, Egam = np.meshgrid(Ex, Egam, indexing='ij')
# 	# Exf = Ex - Egam
# 	# Egam[Egam<50] = 50
# 	# Exf[Exf<300] = 300

# 	# Exf /= 1000   # Convert 
# 	# Egam /= 1000  # to MeV
# 	# Ex /= 1000   # for numerics
# 	# R = np.power(Egam,Nexp)*np.exp(2*np.sqrt(a*Exf))/np.power(Exf,2)
# 	# # print R
# 	# R[np.logical_and(Egam > Ex + (dE_gamma + 0.5*np.abs(Ay1))/1000, Exf > 0)] = 0 # Enforce diagonal requirement same as MAMA
# 	# # print R
# 	# # plt.matshow(R, origin="lower", norm=LogNorm(vmin=0.1, vmax=R.max()))
# 	# # plt.colorbar()
# 	# # plt.show()

# 	# ==== Set up multiplicities ====
# 	multiplicities = np.zeros(2048)
# 	# print
# 	# L1 = np.maximum(np.zeros(Ngates),(slide(Ex,ThresTot,ThresSta,ThresRatio)-Ax0)/Ax1 + 0.5)
# 	# print slide(Ex,ThresTot,ThresSta,ThresRatio)

# 	Ex = Ay0 + Ay1*np.arange(0,Ny,1)
# 	Egamma = Ax0 + Ax1*np.arange(0,Nx,1)

# 	Egamma_weighted_average = np.zeros(Ngates) # The variable called CenS in figega.f
# 	for j in range(Ngates):
# 		j_ex = IyH + j*IyStep # Corresponding Ex channel
# 		Egamma_low_limit = np.maximum(0,slide(Ex[j_ex], ThresTot, ThresSta, ThresRatio))
# 		Egamma_high_limit = Ex[j_ex] + dE_gamma
# 		Egamma_range = np.logical_and(Egamma > Egamma_low_limit, Egamma < Egamma_high_limit)
# 		Egamma_weighted_average[j] = np.sum(Egamma[Egamma_range] * matrix[j_ex, Egamma_range]) / np.sum(matrix[j_ex, Egamma_range])
# 		ExEntry = np.maximum(0, np.minimum(ExEntry0, Ex[j_ex]-200)) # Entry point energy
# 		multiplicities[j] = np.maximum(0, (Ex[j_ex]-ExEntry)/Egamma_weighted_average[j])

# 	plt.plot(multiplicities)
# 	plt.show()
	

# 	# OLD VERSION OF MULTIPLICITY CALCULATION BELOW:
# 	# For each gamma energy gate indexed by j from 0 to Ngates,
# 	# calculate multiplicities:
# 	# (This is adapted from figega.f, I don't understand everything but I 
# 	# think it's identical)
# 	# for j in range(Ngates): # Should try to vectorize away this loop
# 	# 	xm = 0
# 	# 	Exj = Ay0 + Ay1*(IyH + j*IyStep) # Excitation energy
# 	# 	LF2 = (((Exj + dE_gamma - Ax0)/Ax1)+0.5)
# 	# 	if LF2 > len(matrix[0,:]-1):
# 	# 		LF2 = len(matrix[0,:])-1
# 	# 	L1 = np.maximum(0,int((ThresRatio*Exj-Ax0)/Ax1 + 0.5)) #Make sure that slide is used correctly, needs upper threshold!
# 	# 	if StaTot == 2:
# 	# 		L1 = LT1
# 	# 	i_indices = np.arange(L1,LF2,1)
# 	# 	print i_indices
# 	# 	CenS = np.sum(i_indices*matrix[j,i_indices.astype(int)])
# 	# 	summ = np.sum(matrix[j,i_indices.astype(int)])
# 	# 	print CenS
# 	# 	print summ
# 	# 	if summ > 0:
# 	# 		CenS = CenS/summ # Some kind of normalization?
# 	# 	if CenS > 0:
# 	# 		CenS = Ax0 + CenS*Ax1
# 	# 		ExEntry = np.minimum(ExEntry0, Exj-200)
# 	# 		if ExEntry < 0:
# 	# 			ExEntry = 0
# 	# 		xm = (Exj-ExEntry)/CenS
# 	# 	if xm < 0:
# 	# 		xm = 0
# 	# 	Mult[j] = xm

# 	# plt.plot(Mult)
# 	# plt.show()

# 	# Setting up lower limit for weighting function, depends on thresholds
# 	LW1 = int((ThresTot/np.abs(Ay1))+0.5)

# 	# ======================================
# 	# Begin subtraction routine.
# 	# "Fasten seatbelts", as Magne puts it.
# 	# ======================================

# 	# Set up first-iteration weight spectra, assuming normalized squares
# 	R = np.zeros((Ny, Nx))
# 	for i in range(1,Ny):
# 		R[i,0:i] = np.ones(i)/i
# 	Rold = R # This will be used to store previous iteration
# 	# plt.matshow(R, origin='lower')
# 	# plt.colorbar(norm=LogNorm(vmin=0.1, vmax=R.max()))
# 	# plt.show()

# 	fg_matrix = np.copy(matrix) # Make a copy
# 	iterations = 0
# 	while iterations < N_iterations:
# 		for j in range(IyH-IyL): # For each excitation energy, starting from the top, go through procedure to make first-generation spectrum
# 			j_ex = IyH - j # The index of the current excitation energy channel (descending)
# 			Egamma_low_limit = ThresTot # Should this be conditional on StaTot? I.e. ThresSta if Sta?
# 			Egamma_high_limit = Ex[j_ex] + dE_gamma
# 			Egamma_range = np.logical_and(Egamma > Egamma_low_limit, Egamma < Egamma_high_limit)
# 			F = fg_matrix[j_ex,:]
# 			sumF = np.sum(F)
			
# 			# Prepare weight function:
# 			weight = R[j,:]
# 			# Remove negative values and re-normalize
# 			# NOTE TO SELF: Consider implementing more of the features from "Massage"
# 			weight_sum = np.sum(weight)
# 			weight[weight<0] = 0
# 			weight *= weight_sum/np.sum(weight)
# 			# Add a mix of previous iteration to stabilize
# 			if iterations > 4:
# 				weight = 0.7*weight + 0.3*Rold[j,:]
# 			# Again remove negative values and re-normalize
# 			weight_sum = np.sum(weight)
# 			weight[weight<0] = 0
# 			weight *= weight_sum/np.sum(weight)

# 			# Store for next iteration:
# 			Rold[j,:] = weight
			
# 			Gtot = np.zeros(Nx)
# 			for jj in range(j): # Combine all spectra below F to subtract:
# 				jj_ex = j - jj
# 				Egamma_low_limit = np.maximum(0,slide(Ex[jj], ThresTot, ThresSta, ThresRatio)) # Should this be conditional on StaTot? I.e. ThresSta if Sta?
# 				Egamma_high_limit = Ex[jj_ex] + dE_gamma
# 				Egamma_range = np.logical_and(Egamma > Egamma_low_limit, Egamma < Egamma_high_limit)
# 				G = matrix[jj,:]
# 				G[np.logical_not(Egamma_range)] = 0 # Zero everything outside range
# 				sumG = np.sum(G)
# 				denominator = multiplicities[j] * sumG
# 				if denominator > 0:
# 					# Confused about the indices here, should double check!
# 					factor = weight[jj] * multiplicities[jj] * sumF / denominator
# 				else:
# 					factor = 0
# 				Gtot += factor * G

# 			# Calculate area correction alpha:
# 			alpha = 1
# 			sumGtot = np.sum(Gtot)
# 			if sumGtot > 0:
# 				alpha = 1 - 1/multiplicities[j]*sumF/sumGtot
# 			alpha0 = alpha
# 			if alpha0 < 0.85:
# 				alpha0 = 0.85
# 			elif alpha0 > 1.15:
# 				alpha0 = 1.15

# 			# Perform subtraction:
# 			Gtot *= alpha0 # Area correction
# 			FG = F - Gtot
# 			fg_matrix[j_ex,:] = FG

			

# 		R = Rold
		
# 		iterations += 1

# 		plt.matshow(fg_matrix, origin='lower', norm=LogNorm(vmin=0.01, vmax=fg_matrix.max()))
# 		plt.show()





# 	# F = np.zeros(4096)
# 	# SPECi = np.zeros(4096)
# 	# for j in range(IyH, IyL, IyStep):
# 	# 	Ij = np.abs(j-IyH) # Index for Mult and Sing
# 	# 	Exj = Ay0 + Ay1*j  # Excitation energy for uppermost spectrum
# 	# 	LW2 = int(((Exj+dE_gamma-Ay0)/Ay1)+0.5) # Upper limit index for weighting array, taking dE_gamma extra
# 	# 	LW2 = np.abs(LW2 - IyL)
# 	# 	F = matrix[j,:]
# 	# 	print F
# 	# 	# F[]
# 	# 	print F






# # def slide_vectorized(Ex, Thres1, Thres2, ThresRatio):
# # 	slide = ThresRatio*Ex
# # 	# print Thres1, Thres2
# # 	slide[slide<Thres1] = Thres1
# # 	slide[slide>Thres2] = Thres2
# # 	return slide
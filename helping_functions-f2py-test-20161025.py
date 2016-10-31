from pyma_v1 import *

filename = "/home/jorgenem/Dropbox/phd/184W-eksperiment/analysis_187Re/unfolding/alfna-attempted_removing_1p2MeV_Al-20160518.m"
matrix, a, Egamma_range_matrix, Ex_range_matrix = read_mama(filename)
Ex_high = 7500 # keV - maximum excitation energy
Ex_low = 4000
Egamma_low = 1500
Egamma_padding = 900 # We allow gamma energies somewhat higher than Ex_high because of detector uncertainty
dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
N_Exbins = 100
N=120 # Number of points in F and rho

# Get first generation spectrum of exact raw matrix
fgmat, diff_matrix, Egamma_range, Ex_range = first_generation_spectrum_test2(matrix, Egamma_range_matrix, Ex_range_matrix, N_Exbins, Ex_high, dE_gamma, N_iterations=20)
fgvar, a_fgvar, x_fgvar, y_fgvar = read_mama('fgvar_test-20160916-Nstat1000-Nexbins100-corrected_calibration.m')

# Calculate maximal gamma energy we include, given by Ex_high + padding
Egamma_max = Ex_range[-1] + Egamma_padding
if Egamma_max > Egamma_range[-1]:
	Egamma_max = Egamma_range[-1]
# Take out the slice along gamma axis corresponding to these energy limits (also for variance matrix)
fgmat_padded, Egamma_range_sliced = slice_matrix_simple(fgmat, Egamma_range, [0,Egamma_max], axis=1)
fgvar_padded, tmp = slice_matrix_simple(fgvar, Egamma_range, [0,Egamma_max], axis=1)
# plt.figure(0)
# plt.pcolormesh(Egamma_range_sliced, Ex_range, fgmat_padded)
# plt.show()
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

# Normalize the fg matrix and variance matrix (for each Ei bin, which is the way to do it to have normalized branching ratios)
fg_EiEg = div0(fg_EiEg, fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1))
fgv_EiEg = div0(fgv_EiEg, np.power(fg_EiEg.sum(axis=1).reshape(fg_EiEg.shape[0],1),2)) # Normalize the variance accordingly, using that Var(aX) = a^2 Var(X).
# # Need two different versions of matrices, one with Ei and one with Ef on y axis.
fg_EfEg = EitoEf(fg_EiEg, Ex_range_squared)
fgv_EfEg = EitoEf(fgv_EiEg, Ex_range_squared)
# Also get the mask in EfEg coordinates
mask_EfEg = EitoEf(mask_EiEg, Ex_range_squared)

# # Plot masked, squared and rebinned fg to check:
# plt.figure()
# plt.pcolormesh(Eg_range_squared, Ex_range_squared, fg_EiEg)
# plt.show()

	
# Make initial F and rho
rho = np.ones(N)
# F = rebin(fg_EfEg.mean(axis=0), N, rebin_axis=0)
F = fg_EiEg.mean(axis=0) # CHECK THAT THIS IS PROPERLY NORMALIZED





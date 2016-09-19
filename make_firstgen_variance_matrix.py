from pyma_v1 import *

filename = "/home/jorgenem/Dropbox/PhD/184W-eksperiment/analysis_187Re/unfolding/alfna-attempted_removing_1p2MeV_Al-20160518.m"
matrix, a, Egamma_range_matrix, Ex_range_matrix = read_mama(filename)
Ex_max = 7500 # keV - maximum excitation energy
# Ex_min = 300 # keV - minimal excitation energy, effectively moving the ground-state energy up because we cannot resolve the low-energy yrast gamma lines. This is weighed up by also using an effective multiplicity which is lower than the real one, again not considering the low-energy yrast gammas.
dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
# N_Exbins = np.sum(np.logical_and(0 < Ex_range_matrix, Ex_range_matrix < Ex_max + dE_gamma)) # Same number of bins as original matrix
# print "N_Exbins =", N_Exbins
N_Exbins = 100

firstgen_matrix, diff_matrix, Egamma_range_fg, Ex_range_fg = first_generation_spectrum_test2(matrix, Egamma_range_matrix, Ex_range_matrix, N_Exbins, Ex_max, dE_gamma, N_iterations=20)


N_stat = 1000 # How many perturbed copies do we want in our ensemble?
# print np.append(matrix_shape,N_stat)
matrix_ensemble = np.empty(np.append(matrix.shape,N_stat))
firstgen_ensemble = np.empty(np.append(firstgen_matrix.shape,N_stat))

np.random.seed(2)
for i in range(N_stat):
	# matrix_ensemble_current = np.maximum(matrix + np.random.normal(size=matrix_shape)*np.sqrt(matrix), np.zeros(matrix_shape)) # Each bin of the matrix is perturbed with a gaussian centered on the bin count, with standard deviation sqrt(bin count). Also, no negative counts are accepted.
	matrix_ensemble_current = matrix + np.random.normal(size=matrix.shape)*np.sqrt(np.where(matrix > 0, matrix, 0)) # Assuming sigma \approx n^2 / N where n is current bin count and N is total count, according to sigma^2 = np(1-p) for normal approx. to binomial distribution.
	matrix_ensemble_current[matrix_ensemble_current < 0] = 0
	matrix_ensemble[:,:,i] = matrix_ensemble_current
	firstgen_ensemble_current, diff_matrix, tmp2, tmp3 = first_generation_spectrum_test2(matrix_ensemble_current, Egamma_range_matrix, Ex_range_matrix, N_Exbins, Ex_max, dE_gamma, N_iterations=15)
	# if i < 10:
	# 	print "firstgen_ensemble_current.shape =", firstgen_ensemble_current.shape
	firstgen_ensemble[:,:,i] = firstgen_ensemble_current
	# firstgen_matrix, diff_matrix, Egamma_range, Ex_range = first_generation_spectrum_test2(matrix, a, x_array, y_array, N_Exbins, Ex_max, dE_gamma, N_iterations=10)
	print "iteration i =", i

	# Test rhosigchi on current perturbation
	# rhosigchi(firstgen_ensemble_current, fgvar, Egamma_range_fg, Ex_range_fg, N=50, method="Nelder-Mead")

fgvar = firstgen_ensemble.var(axis=2)

# # # Write firstgen matrix and variance to file, to allow for efficient testing
# write_mama(firstgen_matrix, 'firstgen_test.m', Egamma_range_fg, Ex_range_fg)
write_mama(fgvar, 'fgvar_test-20160916-Nstat1000-Nexbins100-corrected_calibration.m', Egamma_range_fg, Ex_range_fg)

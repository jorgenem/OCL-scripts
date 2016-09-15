from pyma_v1 import *

# # Read MAMA matrix
# # filename = "alfna-unfolded-20160518.m"
# filename = "/home/jorgenem/Dropbox/PhD/184W-eksperiment/analysis_187Re/unfolding/alfna-attempted_removing_1p2MeV_Al-20160518.m"
# matrix, a, x_array, y_array = read_mama(filename)
# print a
# Ex_max = 7500 # keV - maximum excitation energy
# # Ex_min = 300 # keV - minimal excitation energy, effectively moving the ground-state energy up because we cannot resolve the low-energy yrast gamma lines. This is weighed up by also using an effective multiplicity which is lower than the real one, again not considering the low-energy yrast gammas.
# dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
# N_Exbins = np.sum(np.logical_and(0 < y_array, y_array < Ex_max + dE_gamma)) # Same number of bins as original matrix
# # N_Exbins = 200

# # #####################################################################################
# # TODO firstgen: 
# # - Why doesn't the convergence criterion work? Does the steady-state oscillate?
# # #####################################################################################

# firstgen_matrix, diff_matrix, Egamma_range, Ex_range = first_generation_spectrum_test2(matrix, a, x_array, y_array, N_Exbins, Ex_max, dE_gamma, N_iterations=10)

# # # print matrix.shape

# # Plot matrix
# plt.figure(0)
# plt.subplot(1,2,1)
# # matrix[matrix < 1] = 0 # I do this to get the colors nice. Should check that it doesn't remove any wanted points, but a fractional count does sound strange (although this is after unfolding)
# plt.pcolormesh(x_array, y_array, matrix, norm=LogNorm(vmin=0.001, vmax=max(matrix.max(), firstgen_matrix.max())), cmap='gist_rainbow_r')
# # plt.matshow(matrix)
# # plt.pcolormesh(matrix, norm=LogNorm(vmin=0.1, vmax=matrix.max()), cmap='gist_rainbow_r')
# # plt.colorbar()
# plt.xlabel('$E_\gamma$ [keV]', fontsize=14)
# plt.ylabel('$E_x [keV]$', fontsize=14)
# # plt.ylim([0,12])
# # plt.xlim([0,12])
# # plt.show()
# # sys.exit(0)
# # END plot matrix



# plt.subplot(1,2,2)
# plt.pcolormesh(Egamma_range, Ex_range, firstgen_matrix, norm=LogNorm(vmin=0.001, vmax=max(matrix.max(), firstgen_matrix.max())), cmap='gist_rainbow_r')
# plt.ylim([0,max(y_array)])
# plt.colorbar()
# plt.xlabel('$E_\gamma$ [keV]', fontsize=14)
# plt.ylabel('$E_x [keV]$', fontsize=14)

# plt.show()


# # # Test writing to file:
# filename_out = "firstgen_test1.m"
# # # write_mama(matrix, filename_out, x_array, y_array)
# write_mama(firstgen_matrix, filename_out, Egamma_range, Ex_range)
# # # print matrix[50:60,50:60]


# sys.exit(0)







##########################################################################
# Test statistical error propagation
##########################################################################
# Read MAMA matrix
# filename = "alfna-unfolded-20160518.m"
filename = "/home/jorgenem/Dropbox/PhD/184W-eksperiment/analysis_187Re/unfolding/alfna-attempted_removing_1p2MeV_Al-20160518.m"
matrix, a, Egamma_range_matrix, Ex_range_matrix = read_mama(filename)
Ex_max = 7500 # keV - maximum excitation energy
# Ex_min = 300 # keV - minimal excitation energy, effectively moving the ground-state energy up because we cannot resolve the low-energy yrast gamma lines. This is weighed up by also using an effective multiplicity which is lower than the real one, again not considering the low-energy yrast gammas.
dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
# N_Exbins = np.sum(np.logical_and(0 < Ex_range_matrix, Ex_range_matrix < Ex_max + dE_gamma)) # Same number of bins as original matrix
# print "N_Exbins =", N_Exbins
N_Exbins = 100

plt.figure(0)
plt.subplot(2,1,1)
plt.pcolormesh(Egamma_range_matrix, Ex_range_matrix, matrix, norm=LogNorm(vmin=0.1, vmax=1e4))

# Get first generation spectrum of exact raw matrix
firstgen_matrix, diff_matrix, Egamma_range_fg, Ex_range_fg = first_generation_spectrum_test2(matrix, Egamma_range_matrix, Ex_range_matrix, N_Exbins, Ex_max, dE_gamma, N_iterations=20)
# firstgen_matrix, a_fg, x_fg, y_fg = read_mama('firstgen_test.m')
# firstgen_matrix, a_fg, x_fg, y_fg = read_mama('firstgen2.tmp')
# print x_fg
# print y_fg
plt.subplot(2,1,2)
plt.pcolormesh(Egamma_range_fg, Ex_range_fg, firstgen_matrix, norm=LogNorm(vmin=0.1, vmax=1e4))
plt.show()

sys.exit(0)

# Get shape of matrices
matrix_shape = matrix.shape
# print "matrix.shape = ", matrix_shape
firstgen_shape = firstgen_matrix.shape
# print "firstgen.shape = ", firstgen_shape


# firstgen2, tmp1, tmp2, tmp3 = first_generation_spectrum_test2(matrix, a, x_array, y_array, N_Exbins, Ex_max, dE_gamma, N_iterations=10)
# print "firstgen2.shape =", firstgen2.shape
# print tmp2
# print tmp3

# write_mama(firstgen2, "firstgen2.m", tmp2, tmp3)

# sys.exit(0)



# Make gaussian random fluctuations
# plt.matshow(matrix, origin='lower')
# plt.colorbar()
N_stat = 1000 # How many perturbed copies do we want in our ensemble?
# print np.append(matrix_shape,N_stat)
matrix_ensemble = np.empty(np.append(matrix_shape,N_stat))
firstgen_ensemble = np.empty(np.append(firstgen_shape,N_stat))

# For rhosigchi testing: import previously made fgvar from file
fgvar, a_fgvar, x_fgvar, y_fgvar = read_mama('fgvar_test-Nstat1000-Nexbins100.m')
# print "fgvar.shape = ", fgvar.shape

# Test rhosigchi using unperturbed FG matrix, importing fgvar from file to use in chisquare denominator:
rhosigchi2(firstgen_matrix, fgvar, Egamma_range_fg, Ex_range_fg, 0, N=100) # rhosigchi2 is a code-up of Schiller's original method
# rhosigchi(firstgen_matrix, fgvar, Egamma_range_fg, Ex_range_fg, N=50, method="Nelder-Mead")
sys.exit(0)


# # TODO: Parallellize this. Or maybe it's not needed? Seems Numpy auto-parallelizes the hard tasks :D
np.random.seed(2)
for i in range(N_stat):
	# matrix_ensemble_current = np.maximum(matrix + np.random.normal(size=matrix_shape)*np.sqrt(matrix), np.zeros(matrix_shape)) # Each bin of the matrix is perturbed with a gaussian centered on the bin count, with standard deviation sqrt(bin count). Also, no negative counts are accepted.
	matrix_ensemble_current = matrix + np.random.normal(size=matrix_shape)*np.sqrt(np.where(matrix > 0, matrix, 0)) # Assuming sigma \approx n^2 / N where n is current bin count and N is total count, according to sigma^2 = np(1-p) for normal approx. to binomial distribution.
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
# print "Freshly calculated fgvar.shape =", fgvar.shape
# plt.matshow(fgvar, origin='lower')
# plt.colorbar()


# plt.show()



# # # Write firstgen matrix and variance to file, to allow for efficient testing
# write_mama(firstgen_matrix, 'firstgen_test.m', Egamma_range_fg, Ex_range_fg)
# write_mama(fgvar, 'fgvar_test-Nstat1000-Nexbins100.m', Egamma_range_fg, Ex_range_fg)





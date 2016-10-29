from pyma_v1 import *

filename = "/home/jorgenem/Dropbox/phd/184W-eksperiment/analysis_187Re/unfolding/alfna-attempted_removing_1p2MeV_Al-20160518.m"
matrix, a, Egamma_range_matrix, Ex_range_matrix = read_mama(filename)
Ex_high = 7500 # keV - maximum excitation energy
Ex_low = 4000
Egamma_low = 1500
Egamma_padding = 900 # We allow gamma energies somewhat higher than Ex_high because of detector uncertainty
dE_gamma = 300 # keV - allow gamma energy to exceed excitation energy by this much, to account for experimental resolution
N_Exbins = 100
N=150 # Number of points in F and rho

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

# Make Ex and Eg range arrays including padding
# Note 20161029: I added this to be able to use rebin_and_shift() instead of just rebin(), and thus get calibrations right. 
# However, I realize that since rebin_and_shift does not support below-zero energy range, it will not work just like this.
# I think it can be done by reordering the above, and first rebinning on the non-padded (possibly only padded to positive Ex) 
# array and then pad (it should be even easier this way, since calibration is then automatically the same for both axes).
Ex_range_padded = np.linspace(Ex_range[0]-N_Ex_padding_bottom*(Ex_range[1]-Ex_range[0]), Ex_range[-1]+N_Ex_padding_top*(Ex_range[1]-Ex_range[0]), fgmat_padded.shape[0])
Eg_range_padded = np.linspace(Egamma_range_sliced[0]-N_Eg_padding_bottom*(Egamma_range_sliced[1]-Egamma_range_sliced[0]), Egamma_range_sliced[-1], fgmat_padded.shape[1])

print "Ex_range_padded =",Ex_range_padded
print "Eg_range_padded =",Eg_range_padded
sys.exit(0)

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

for iteration in range(20):
	if iteration < 4:
		deviation_limit = 1.2
	elif iteration < 11:
		deviation_limit = 1.1
	elif iteration < 20:
		deviation_limit = 1.05
	elif iteration < 29:
		deviation_limit = 1.025
	else:
		deviation_limit = 1.01	


	# === Calculate helping functions s, a and b: ===
	
	# Allocate arrays to be filled
	s = np.zeros(N)
	a = np.zeros(N)
	b = np.zeros(N)
	
	# Find the index corresponding to lower Eg and lower/upper Ei limits:
	i_Eg_min = np.where(Eg_range_squared > Egamma_low)[0][0]
	if np.abs(Eg_range_squared[i_Eg_min]) > np.abs(Eg_range_squared[i_Eg_min-1]):
		i_Eg_min -= 1
	i_Ei_min = np.where(Ex_range_squared > Ex_low)[0][0]
	if np.abs(Ex_range_squared[i_Ei_min]) > np.abs(Ex_range_squared[i_Ei_min-1]):
		i_Ei_min -= 1
	i_Ei_max = np.where(Ex_range_squared > Ex_high)[0][0]
	if np.abs(Ex_range_squared[i_Ei_max]) > np.abs(Ex_range_squared[i_Ei_max-1]):
		i_Ei_max -= 1
	
	# Fill the arrays:
	for i_Ei in range(i_Ei_min,i_Ei_max):
		s[i_Ei] = ( F[i_Eg_min:i_Ei] * np.flipud(rho[0:i_Ei-i_Eg_min]) ).sum()
		b[i_Ei] = ( div0( F[i_Eg_min:i_Ei] * np.flipud(rho[0:i_Ei-i_Eg_min]) * fg_EiEg[i_Ei,i_Eg_min:i_Ei], np.power(fgv_EiEg[i_Ei,i_Eg_min:i_Ei],2) ) ).sum()
		a[i_Ei] = ( np.power( div0( F[i_Eg_min:i_Ei] * np.flipud(rho[0:i_Ei-i_Eg_min]), fgv_EiEg[i_Ei,i_Eg_min:i_Ei]), 2) ).sum()
	

	# plt.figure()
	# plt.title('s,a,b')
	# plt.plot(s, label='s')
	# plt.plot(b, label='b')
	# plt.plot(a, label='a')
	# plt.yscale('log')
	# plt.legend()
	
	# === Make psi and phi matrices: ===
	
	# # Smooth the fgv matrix
	# from scipy.signal import savgol_filter
	# fgv_EiEg = savgol_filter(fgv_EiEg, 5, 3, axis=1)
	# # Apply mask to smoothed fgv
	# fgv_EiEg = mask_EiEg*fgv_EiEg
	
	# Remove very low values of fgv
	# fgv_EiEg[fgv_EiEg<1e-2] = 0
	
	psi = div0( 1, np.power(s.repeat(N).reshape((N,N)) * fgv_EiEg, 2) ) # Plotting psi makes it seem that there are just a few non-zero values, but using matshow and inspecting the values I see that in fact all values are resolved -- they are just very small. But absent numerical instabilities, this is what the recipe prescribes, so I will code it up this way and see what happens.
	phi = (div0(a,np.power(s,3))).repeat(N).reshape(N,N) - (div0(b,np.power(s,2))).repeat(N).reshape((N,N)) + div0(fg_EiEg, s.repeat(N).reshape((N,N))*np.power(fgv_EiEg, 2))
	phi = phi*mask_EiEg # This one needs masking
	
	# plt.matshow(phi, origin='lower')
	# plt.colorbar()
	# plt.show()
	
	
	# === Calculate F and rho: ===
	
	F_new = np.zeros(N)
	for i_Eg in range(i_Eg_min, i_Ei_max):
		F_new_current = (rho[max(i_Ei_min-i_Eg,0):i_Ei_max-i_Eg] * phi[max(i_Ei_min,i_Eg):i_Ei_max,i_Eg]).sum() / (np.power(rho[max(i_Ei_min-i_Eg,0):i_Ei_max-i_Eg],2) * psi[max(i_Ei_min,i_Eg):i_Ei_max,i_Eg]).sum()
		if F_new_current > deviation_limit*F[i_Eg]:
			F_new[i_Eg] = deviation_limit*F[i_Eg]
		elif	F_new_current < F[i_Eg]/deviation_limit:
			F_new[i_Eg] = F[i_Eg]/deviation_limit
		else:
			F_new[i_Eg] = F_new_current



	# print np.isnan(F_new)
	# Since we're not using div0() function for F_new, we need to ensure no NaNs were produced.
	F_new[np.isnan(F_new)] = 0
	# print np.isnan(F_new)
	

	rho_new = div0( (F.repeat(N).reshape((N,N)).T * EitoEf(phi, Ex_range_squared)).sum(axis=1), ( (np.power(F,2)).repeat(N).reshape((N,N)).T * EitoEf(psi, Ex_range_squared) ).sum(axis=1) )
	rho_new[rho_new > deviation_limit*rho] = deviation_limit*rho[rho_new > deviation_limit*rho]
	rho_new[rho_new < deviation_limit*rho] = rho[rho_new < deviation_limit*rho]/deviation_limit
	# A more explicit calculation of rho. UPDATE 20161028: Didn't finish this, going back to trying the above version.
	# rho_new = np.zeros(N)
	# for i_Ef in range(0, i_Ei_max):
	# 	rho_new[i_Ef] = F[]

	
	plt.figure()
	plt.title('F_new / rho_new')
	plt.plot(F_new, label='F_new')
	plt.plot(rho_new, label='rho_new')
	plt.plot(F, label='F_old')
	plt.plot(rho, label='rho_old')
	plt.legend()
	plt.yscale('log')
	# plt.show()
	
	# Calculate fgtheo
	plt.figure()
	fgtheo = F_new.repeat(N).reshape((N,N)).T * EftoEi(rho_new.repeat(N).reshape((N,N)), Ex_range_squared) * mask_EiEg
	# fgtheo_denom = fgtheo.sum(axis=1).repeat(N).reshape((N,N)).T
	fgtheo_denom = s.repeat(N).reshape((N,N)).T
	fgtheo = div0(fgtheo, fgtheo_denom)
	plt.pcolormesh(fgtheo)
	plt.colorbar()
	plt.title('fgtheo')
	plt.show()

	# Switch naming for next iteration
	# F_old = F
	# rho_old = rho
	F = F_new.copy()
	rho = rho_new.copy()
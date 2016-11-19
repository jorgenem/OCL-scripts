import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/OCL-scripts')
import pyma_v1 as pyma
# sys.path.insert(0,'/home/jorgenem/gitrepos/oslo-method-software/prog')
import rhosigchi_f2py as rsc
import numpy as np 
import matplotlib.pyplot as plt

# print rsc.__doc__
# sys.exit(0)

# Setting limits on gamma and excitation energy
Eg_min = 960
Ex_min = 4000
Ex_max = 7000

# Read main fg and fgv:
fg, calib, Eg_range, Ex_range = pyma.read_mama('firstgen-jem-20161011-Nexbins196-Nstat500.m')
fgv, tmp1, tmp2, tmp3 = pyma.read_mama('fgvar-jem-20161011-Nexbins196-Nstat500.m')

# Reformat: Put fg and fgv into 512x512 matrices to comply with rhosigchi:
fg_reformatted = np.zeros((512,512))
fg_reformatted[0:fg.shape[0],0:fg.shape[1]] = fg
fgv_reformatted = np.zeros((512,512))
fgv_reformatted[0:fgv.shape[0],0:fgv.shape[1]] = fgv
fg_reformatted = fg_reformatted.T 	# Fortran uses 
fgv_reformatted = fgv_reformatted.T # column-major order

# Calculate main rho and T
rho, T, calib_out = rsc.rhosigchi(fg_reformatted, fgv_reformatted, calib, Eg_min, Ex_min, Ex_max)
print calib_out
sys.exit(0)

#### CHOICE 1: Make rsc ensemble 
# Read in ensemble of perturbed fg matrices and calculate rho and T for each:
N_ensemble = 1
rho_var_array = np.zeros((N_ensemble,101))
T_var_array = np.zeros((N_ensemble,101))
for i in range(N_ensemble):
	filename = "fg_ensemble/firstgen_ensemble-Nexbins196-%03d.m"%i
	fg_current, tmp1, tmp2, tmp3 = pyma.read_mama(filename)

	# Reformat: 
	fg_reformatted = np.zeros((512,512))
	fg_reformatted[0:fg.shape[0],0:fg.shape[1]] = fg_current
	fg_reformatted = fg_reformatted.T 	
	# Calculate current rho and T
	rho_var_array[i,:], T_var_array[i,:] = rsc.rhosigchi(fg_reformatted, fgv_reformatted, calib, Eg_min, Ex_min, Ex_max)

sys.exit(0)

# # Save ensemble arrays to file
# np.savetxt("rsc_rho_var_array-Nensemble%03d"%N_ensemble.dat, rho_var_array)
# np.savetxt("rsc_T_var_array-Nensemble%03d"%N_ensemble.dat, T_var_array)

#### CHOICE 2: Read previously generated rsc ensemble from file
rho_var_array = np.loadtxt("rsc_rho_var_array-Nensemble500.dat")
T_var_array   = np.loadtxt("rsc_T_var_array-Nensemble500.dat")
N_ensemble = rho_var_array.shape[0]

# == Plotting ==
# Set up axis arrays
E_rho = np.linspace(0,100,101)*120 - 840
E_T   = np.linspace(0,100,101)*120 + 960

# plt.figure(1)
# plt.title('Ensemble perturbations: Rho')
# plt.plot(E_rho, rho)
# for i in range(N_ensemble):
# 	plt.plot(E_rho, rho_var_array[i,:])
# plt.yscale('log')

# plt.figure(2)
# plt.title('Ensemble perturbations: T')
# plt.plot(E_T, T)
# for i in range(N_ensemble):
# 	plt.plot(E_T, T_var_array[i,:])
# plt.yscale('log')
# plt.show()



# === Error bar/band calculations ===
plt.figure(5)
Nbins = 50
rho_hist_points_list = [20,40,50]
for i in rho_hist_points_list:
    rhohistcurrent, tmp = np.histogram(rho_var_array[:,i], bins=Nbins)
    plt.step(np.linspace(0,Nbins,Nbins), rhohistcurrent, label='$E_x$ = {}'.format(E_rho[i]), linewidth=1.5)
plt.legend()    
plt.savefig('rho_hist.png')


plt.figure(6)
Nbins = 50
T_hist_points_list = [20,40,60]
for i in T_hist_points_list:
    Thistcurrent, tmp = np.histogram(T_var_array[:,i], bins=Nbins)
    plt.step(np.linspace(0,Nbins,Nbins), Thistcurrent, label='$E_\gamma$ = {}'.format(E_T[i]), linewidth=1.5)
plt.legend()    

plt.savefig('T_hist.png')



# sys.exit(0)



percentiles = [2.5,5,15.87,50,85.14,95,97.5]

plt.figure(1)
rhopercentiles = np.percentile(rho_var_array, percentiles, axis=0)
# for i in range(rhopercentiles.shape[0]):
# 	plt.plot(E_rho, rhopercentiles[i,:], label="percentile %d"%percentiles[i])
plt.plot(E_rho, rho, label='original', color='k', linewidth=1)
plt.fill_between(E_rho, rhopercentiles[0,:], rhopercentiles[6,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='blue',
    linewidth=0,label='95\% confidence')
plt.fill_between(E_rho, rhopercentiles[1,:], rhopercentiles[5,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='red',
    linewidth=0,label='90\% confidence')
plt.fill_between(E_rho, rhopercentiles[2,:], rhopercentiles[4,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='gray',
    linewidth=0,label='69\% confidence')
plt.title('Rho percentiles')
plt.xlabel('$E_x$')
plt.yscale('log')
plt.legend()
plt.savefig('rho_with_error_bands.png')


plt.figure(2)
Tpointdist = T_var_array[:,50]
# plt.hist(Tpointdist)
Tpercentiles = np.percentile(T_var_array, percentiles, axis=0)
# for i in range(Tpercentiles.shape[0]):
# 	plt.plot(E_T, Tpercentiles[i,:], label="percentile %d"%percentiles[i])
plt.plot(E_T, T, label='original', color='k', linewidth=1)
plt.fill_between(E_T, Tpercentiles[0,:], Tpercentiles[6,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='blue',
    linewidth=0,label='95\% confidence')
plt.fill_between(E_T, Tpercentiles[1,:], Tpercentiles[5,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='red',
    linewidth=0,label='90\% confidence')
plt.fill_between(E_T, Tpercentiles[2,:], Tpercentiles[4,:],
    alpha=0.5, edgecolor='#3F7F4C', facecolor='gray',
    linewidth=0,label='69\% confidence')
plt.title('T percentiles')
plt.xlabel('$E_\gamma$')
plt.yscale('log')
plt.legend()
plt.savefig('T_with_error_bands.png')



plt.figure(3)
Nbins = 40
rhohist = np.zeros((Nbins,101))
for i in range(101):
    rhohist_curr, bins_curr = np.histogram(rho_var_array[:,i], bins=Nbins)
    rhohist[:,i] = rhohist_curr

plt.pcolormesh(E_rho[5:57],np.linspace(0,Nbins,Nbins),rhohist[:,5:57])
plt.xlabel('$E_x$')
plt.title('Rho: Distribution of ensemble values for each $E_x$')
plt.savefig('rho_2D_histogram.png')


plt.figure(4)
Nbins = 40
Thist = np.zeros((Nbins,101))
for i in range(101):
    Thist_curr, bins_curr = np.histogram(T_var_array[:,i], bins=Nbins)
    Thist[:,i] = Thist_curr

plt.pcolormesh(E_T[15:67],np.linspace(0,Nbins,Nbins),Thist[:,15:67])
plt.xlabel('$E_\gamma$')
plt.title('T: Distribution of ensemble values for each $E_\gamma$')
plt.savefig('T_2D_histogram.png')

# plt.show()
from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import sys

# Constants for energy binning
a0 =  -0.8400
a1 =   0.1200

# Cross section factor for formula strength = xsec*factor/Egamma
xsec_factor = 8.68e-8



# Import data
fermigasfile = open('fermigas.cnt')
fermigaslines = fermigasfile.readlines()
fermigasfile.close()
Nfermi = len(fermigaslines)
energy = np.zeros(Nfermi)
fermigas = np.zeros(Nfermi)
for i in range(Nfermi):
	energy[i] = a0 + i*a1
	fermigas[i] = float(fermigaslines[i].split()[0])

rholevfile = open('rholev.cnt')
rholevlines = rholevfile.readlines()
rholevfile.close()
Nrholev = len(rholevlines)
rholev = np.zeros(Nrholev)
for i in range(Nrholev):
	rholev[i] = float(rholevlines[i].split()[0])

rholevfile = open('rholev.cnt')
rholevlines = rholevfile.readlines()
rholevfile.close()
Nrholev = len(rholevlines)
rholev = np.zeros(Nrholev)
for i in range(Nrholev):
	rholev[i] = float(rholevlines[i].split()[0])

rhopawfile = open('rhopaw.cnt')
rhopawlines = rhopawfile.readlines()
rhopawfile.close()
Nrhopaw = len(rhopawlines)
rhopaw = np.zeros((Nrhopaw,2))
for i in range(Nrhopaw):
	if i < int(Nrhopaw/2):
		rhopaw[i,0] = float(rhopawlines[i].split()[0])
	else:
		rhopaw[i-int(Nrhopaw/2),1] = float(rhopawlines[i].split()[0])

# Level density at neutron separation energy
Bn=7.359200
Bnerr=0.001
rho_Bn=14875000.000000
rho_Bnerr=363000.000000

# Plotting, level density
plt.figure(10)
plt.plot(energy[energy>0], fermigas[energy>0], '--', color='grey', label='Constant-temperature model')
plt.hold('on')
plt.plot(energy[0:Nrholev], rholev, color='black', label='Known levels')
plt.errorbar(energy[0:Nrhopaw], rhopaw[:,0], yerr=rhopaw[:,1], fmt='.', color='midnightblue', label='Present work, exp. data points')
plt.plot(Bn, rho_Bn, 's', label=r'$\rho$ from neutron res. data', color='black')
plt.yscale('log')
plt.xlim([-1,8])
plt.ylim([1e-0, 1e8])
plt.legend(loc='upper left', fontsize=13)
plt.ylabel(r'Level density $\rho (E)$ [MeV$^{-1}$]', fontsize=13)
plt.xlabel('Excitation energy [MeV]', fontsize=13)
plt.text(0, 1e5, '$^{187}\mathrm{Re}$', fontsize=30)
plt.text(-1,1e6, 'PRELIMINARY', alpha=0.1, fontsize=70, rotation=30)
plt.savefig('level_density_pyplot.png')
# plt.show()
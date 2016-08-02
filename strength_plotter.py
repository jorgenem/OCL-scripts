from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import sys

# Constants for energy binning
a0 =  -0.8400
a1 =   0.1200

# Cross section factor for formula strength = xsec*factor/Egamma
xsec_factor = 8.68e-8

strengthfile = open('strength.nrm', 'r')
strengthlines = strengthfile.readlines()
strengthfile.close()
N = len(strengthlines)
strength = np.zeros((N,3))
for i in range(N):
	words = strengthlines[i].split()
	if i < int(N/2):
		strength[i,0] = a0 + i*a1 # Energy coordinate 
		strength[i,1] = float(words[0]) # Strength coordinate
	else:
		strength[i-int(N/2),2] = float(words[0]) # Strength uncertainty (this way due to horrible file format!)

transextfile = open('transext.nrm')
transextlines = transextfile.readlines()
transextfile.close()
Next = len(transextlines)
strengthext = np.zeros((Next, 2))
for i in range(Next):
	transext_current = float(transextlines[i].split()[0])
	energy_current = a0 + i*a1 + 1e-8
	strengthext[i,:] = ( energy_current, transext_current/(2*np.pi*energy_current**3) )

comptonfile = open('187Re_gamma_n.dat', 'r')
comptonlines = comptonfile.readlines()
comptonfile.close()
Ncompton = len(comptonlines)
comptonstrength = np.zeros((Ncompton,3))
for i in range(Ncompton):
	words = comptonlines[i].split()
	comptonstrength[i,0] = float(words[3])
	comptonstrength[i,1] = float(words[0]) * xsec_factor / comptonstrength[i,0]
	comptonstrength[i,2] = np.sqrt(float(words[1])**2 + float(words[2])**2) * xsec_factor / comptonstrength[i,0] # Energy, strength, uncertainty (sqrt(stat^2 + sys^2))



# Plotting:
plt.figure(0)
plt.plot(strengthext[10:,0], strengthext[10:,1], color='navy', label='Present work, extrapolated')
plt.errorbar(strength[:,0], strength[:,1], yerr=strength[:,2], label="Present work, exp. data points", fmt='.', color='navy')
plt.errorbar(comptonstrength[:,0], comptonstrength[:,1], yerr=comptonstrength[:,2], label='Shizuma et.al. (2005)', fmt='.-', color='crimson')
plt.yscale('log')
plt.legend(loc='upper left', fontsize=13)
plt.xlabel('$\gamma$-ray energy $E_\gamma$ [MeV]', fontsize=13)
plt.ylabel('$\gamma$-ray strength [MeV$^{-3}$]', fontsize=13)
plt.text(0, 1e-7, '$^{187}\mathrm{Re}$', fontsize=30)
plt.text(-1.5,1e-6, 'PRELIMINARY', alpha=0.1, fontsize=70, rotation=30)
plt.savefig('strength_pyplot.png')




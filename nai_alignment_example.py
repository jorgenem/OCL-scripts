from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import os
from OCL_calibration_library import *
import ROOT as R
import sys
# Example file for Python function library for lining up and calibrating OCL (and other) detectors
# Written by J{\o}rgen E. Midtb{\o}
# Spring 2016
# j.e.midtbo@fys.uio.no (feel free to send me questions!)
# http://folk.uio.no/jorgeem
# github.com/jorgenem
#

# Read NaI energy histograms, directly from a ROOT file. The name of the ROOT spectrum is m_nai_e (which can be changed).
file1 = R.TFile.Open("example_spectra-NaI.root")
nai = []
nai_e_hist0 = file1.m_nai_e.ProjectionX("myHist%d" %0, 1,1)
xmin = nai_e_hist0.GetXaxis().GetXmin()
xmax = nai_e_hist0.GetXaxis().GetXmax()
Nbins = nai_e_hist0.GetSize()
x_range = np.linspace(xmin, xmax, Nbins+1)
for j in range(32):
    nai_current = np.zeros((4, Nbins))
    nai_current[0,:] = np.linspace(1,Nbins,Nbins)
    nai_current[1,:] = x_range[0:-1]
    nai_current[2,:] = x_range[1:]
    nai_current[3,:] = np.array(file1.m_nai_e.ProjectionX("myHist%d" %j, j+1,j+1))
    nai.append(nai_current)

    
#print nai[0]


# # Plot to check
#%matplotlib auto
# detector = 0
# plt.figure(10)
# plt.hold('on')
# for i in range(27):
# 	strip_current = nai[i]
# 	strip_current = rebin_data(strip_current, 4)
# 	plt.step((strip_current[1,:]+strip_current[2,:])/2,strip_current[3,:], label="spectrum %d" %i)
# plt.legend(fontsize="small")
# plt.show()



# Align individual detectors
i_current = 3
show_error_plots = True
#		0			1			2			3			4			5			6			7			8			9			10			11			12			13			14			15			16			17			18			19			20			21			22			23			24			25			26			27			28			29			30			31
indices=[1, 		1, 			1, 			1, 			0, 			1, 			1, 			1, 			0, 			0, 			1, 			1, 			0, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			1, 			0,			0] # Was detector on? True/False
Nbins = [[0], 		[500],		[500],		[500],		[0], 		[500],		[500],		[500],		[0],		[700],		[500],		[450],		[0],		[450],		[450],		[450],		[450],		[450],		[450],		[450],		[450],		[100,250],	[450],		[350],		[550],		[450],		[450],		[550],		[275],		[275],		[0],		[0]]
binshift=[[0],		[100], 		[110],		[110],		[0], 		[110], 		[110], 		[110], 		[0], 		[40], 		[110], 		[110], 		[0], 		[110], 		[110], 		[110], 		[110], 		[110],		[110], 		[110], 		[110], 		[110,400],	[110], 		[55], 		[110], 		[100], 		[110], 		[110], 		[55], 		[55], 		0,			0]
rebin = [1, 		2, 			2, 			2, 			1, 			2, 			2, 			2, 			0, 			1, 			2, 			2, 			1, 			2, 			2, 			2, 			2, 			2, 			2, 			2, 			2, 			2, 			2, 			4, 			2, 			2, 			2, 			2, 			4, 			4, 			0,			0]
shiftw=[50,			50, 		50, 		50, 		0, 			50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		20, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		50, 		0,			0]
erropt =[1,			3,			3,			3,			1,			3,			3,			3,			0,			1,			3,			3,			1,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			3,			0,			0]
line_up_nai(nai, i_current, i_reference=0, Nbins=Nbins[i_current], binshift=binshift[i_current], rebin=rebin[i_current], shift_width=shiftw[i_current], error_option=erropt[i_current], normalize_spectra=True, show_plots=True, show_error_plots=show_error_plots)




# =========================================================================
sys.exit(0) # Comment out this to run through all detectors at once.
# =========================================================================





gains_list, shifts_list = calibrate_nai_multiple_regions(nai, indices, Nbins, rebin, binshift, shift_width=shiftw, error_option=erropt, print_on=False)

# Print gains and shifts nice
print "Gains:"
for i in range(4):
	i *= 7
	print "%f %f %f %f %f %f %f" %(gains_list[i+0],gains_list[i+1],gains_list[i+2],gains_list[i+3],gains_list[i+4],gains_list[i+5],gains_list[i+6]) 		
print "%f %f %f %f" %(gains_list[28],gains_list[29],gains_list[30],gains_list[31])
print "Shifts:"
for i in range(4):																
	i *= 7
	print "%f %f %f %f %f %f %f" %(shifts_list[i+0],shifts_list[i+1],shifts_list[i+2],shifts_list[i+3],shifts_list[i+4],shifts_list[i+5],shifts_list[i+6]) 		
print "%f %f %f %f" %(shifts_list[28],shifts_list[29],shifts_list[30],shifts_list[31])
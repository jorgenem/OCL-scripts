from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import os
from OCL_calibration_library import *
import ROOT as R
import sys
# Example file for Python function library for lining up and calibrating OCL (and other) detector spectra
# Written by J{\o}rgen E. Midtb{\o}
# Spring 2016
# j.e.midtbo@fys.uio.no (feel free to send me questions!)
# http://folk.uio.no/jorgeem
# github.com/jorgenem
#


# For E-dE alignment and calibration it is fruitful to apply 2D cuts to the spectra before 
# projecting to 1D histograms. At present this has to be done elsewhere - typically in the
# user_sort.cpp routine where the .root file is made. The best option is then to use ROOT's
# TCut objects and methods to apply the cuts. See my example user_sort.cpp file in the Github
# repo for an example of how to implement this.



# Read and project E and dE histograms, directly from a ROOT file.
# This process is quite CPU-heavy, so it should either be saved to a text file or this script
# should be run in an interactive environment like iPython (Notebook) to store them in memory.
# file1 = R.TFile.Open("example_spectra-EdE.root") # Example without TCut
file1 = R.TFile.Open("example_spectra-EdE-TCut.root") # Example with TCut
front = []
back = []
# Get x axis limits
front_b0f0 = file1.m_e_de_b0f0.ProjectionY("myHist_f%d" %0)
xmin_front = front_b0f0.GetXaxis().GetXmin()
xmax_front = front_b0f0.GetXaxis().GetXmax()
Nbins_front = front_b0f0.GetSize()
x_range_front = np.linspace(xmin_front, xmax_front, Nbins_front+1)
back_b0f0 = file1.m_e_de_b0f0.ProjectionX("myHist_b%d" %0)
xmin_back = back_b0f0.GetXaxis().GetXmin()
xmax_back = back_b0f0.GetXaxis().GetXmax()
Nbins_back = back_b0f0.GetSize()
x_range_back = np.linspace(xmin_back, xmax_back, Nbins_back+1)
for i in range(8):
	for j in range(8):
		front_current = np.zeros((4, Nbins_front))
		front_current[0,:] = np.linspace(1,Nbins_front,Nbins_front)
		front_current[1,:] = x_range_front[0:-1]
		front_current[2,:] = x_range_front[1:]
		front_current[3,:] = np.array(file1.Get("m_e_de_b%df%d" %(i,j)).ProjectionY("myHist_f%d" %j))
		front.append(front_current)

		back_current = np.zeros((4, Nbins_back))
		back_current[0,:] = np.linspace(1,Nbins_back,Nbins_back)
		back_current[1,:] = x_range_back[0:-1]
		back_current[2,:] = x_range_back[1:]
		back_current[3,:] = np.array(file1.Get("m_e_de_b%df%d" %(i,j)).ProjectionX("myHist_b%d" %i))
		back.append(back_current)



	
# # Plot to check
# plt.figure(0)
# plt.hold('on')
# strip = 0
# for i in range(8):
# 	strip_current = front[8*i+strip]
# 	strip_current = rebin_data(strip_current, 4)
# 	plt.step((strip_current[1,:]+strip_current[2,:])/2,strip_current[3,:])
# plt.show()

# # ========================================
# sys.exit(0) # Comment out to proceed down
# # ========================================



# Calibrate a single strip of front/back to determine which binshift (and possibly Nbins) to use.
Nbins = [[[35, 100],[35, 100],[35, 100],[35, 100],[35, 100],[35, 100],[35, 100],[35, 100]], 
		 [[150],[150],[150],[150],[150],[150],[150],[150]]]
binshift=[[[70, 122],[70, 135],[70, 140],[70, 140],[70, 140],[70, 140],[70, 140],[70, 135]],
		  [[0],[0],[0],[0],[0],[0],[0],[0]]]
rebin = [[4,4,4,4,4,4,4,4],
		 [4,4,4,4,4,4,4,4]]
# error_option = [[1,1,1,1,1,1,1,1],
# 				[1,1,1,1,1,1,1,1]]
error_option = [[3,3,3,3,3,3,3,3],
				[3,3,3,3,3,3,3,3]]
shift_width = [	[20, 20, 20, 20, 20, 20, 20, 20],
				[20, 20, 20, 20, 20, 20, 20, 20]]



k = 0
if k==0:
	spectra=back
elif k==1:
    spectra=front
strip = 7
g,s = calibrate_ede_strip(spectra, strip, Nbins=Nbins[k][strip], binshift=binshift[k][strip], rebin=rebin[k][strip], error_option=error_option[k][strip], shift_width=shift_width[k][strip], show_plots=True, normalize_spectra=True)
print g,s


# ========================================
sys.exit(0) # Comment out to proceed down
# ========================================


# After finding optimal parameters for each strip, run through the whole thing
back_front_counter = 0
for spectra in [back, front]:
	if back_front_counter == 0:
		print "Back:"
		k = 0
	else:
		print "Front:"
		k = 1
	gains_list = []
	shifts_list = []
	for strip in range(8):
		gains, shifts = calibrate_ede_strip(spectra, strip, Nbins=Nbins[k][strip], rebin=rebin[k][strip], error_option=error_option[k][strip], binshift=binshift[k][strip], shift_width=shift_width[k][strip], normalize_spectra=True, show_plots=True)
		gains_list.append(gains)
		shifts_list.append(shifts)
	print "Gains:"
	for i in range(8):
		# All entries on the same row should be the same physical detector, columns are the same ring/strip/angle
		# The gains and shifts are on the other hand found for each strip, which means all different physical detectors																		
		print "%f %f %f %f %f %f %f %f" %(gains_list[0][i], gains_list[1][i], gains_list[2][i], gains_list[3][i], gains_list[4][i], gains_list[5][i], gains_list[6][i], gains_list[7][i]) 		
	print "Shifts:"
	for i in range(8):
		# All entries on the same row should be the same physical detector, columns are the same ring/strip/angle
		# The gains and shifts are on the other hand found for each strip, which means all different physical detectors																		
		print "%f %f %f %f %f %f %f %f" %(shifts_list[0][i], shifts_list[1][i], shifts_list[2][i], shifts_list[3][i], shifts_list[4][i], shifts_list[5][i], shifts_list[6][i], shifts_list[7][i]) 		
	back_front_counter += 1

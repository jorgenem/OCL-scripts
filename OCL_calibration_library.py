from __future__ import division
# Python function library for lining up and calibrating OCL (and other) detectors
# Written by J{\o}rgen E. Midtb{\o}
# Spring 2016
# j.e.midtbo@fys.uio.no (feel free to send me questions!)
# http://folk.uio.no/jorgeem
# github.com/jorgenem
#

import numpy as np
import sys
from matplotlib import pyplot as plt
from scipy.signal import find_peaks_cwt
from scipy.optimize import minimize
import ROOT as R




def read_and_project_ede_spectra(root_filename):
	file1 = R.TFile.Open(root_filename) 
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
			# if i==0:
			# 	back_current[3,:] /= 10
			back.append(back_current)
	return front, back

def read_and_project_ede_spectra_subregion(root_filename, front_bins, back_bins):
	# Assumes front_bins and back_bins are lists of first and last bin number.
	file1 = R.TFile.Open(root_filename) 
	front = []
	back = []
	# Get x axis limits
	front_b0f0 = file1.m_e_de_b0f0.ProjectionY("myHist_f%d" %0, front_bins[0], front_bins[1])
	xmin_front = front_b0f0.GetXaxis().GetXmin()
	xmax_front = front_b0f0.GetXaxis().GetXmax()
	Nbins_front = front_b0f0.GetSize()
	x_range_front = np.linspace(xmin_front, xmax_front, Nbins_front+1)
	back_b0f0 = file1.m_e_de_b0f0.ProjectionX("myHist_b%d" %0, back_bins[0], back_bins[1])
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
			front_current[3,:] = np.array(file1.Get("m_e_de_b%df%d" %(i,j)).ProjectionY("myHist_f%d" %j, front_bins[0], front_bins[1]))
			front.append(front_current)
	
			back_current = np.zeros((4, Nbins_back))
			back_current[0,:] = np.linspace(1,Nbins_back,Nbins_back)
			back_current[1,:] = x_range_back[0:-1]
			back_current[2,:] = x_range_back[1:]
			back_current[3,:] = np.array(file1.Get("m_e_de_b%df%d" %(i,j)).ProjectionX("myHist_b%d" %i, back_bins[0], back_bins[1]))
			back.append(back_current)
	return front, back

def read_nai_spectra(root_filename):
	file1 = R.TFile.Open(root_filename)
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
	return nai

def read_histogram(root_filename, histogram_object_name):
	file1 = R.TFile.Open(root_filename)
	xmin = file1.Get(histogram_object_name).GetXaxis().GetXmin()
	xmax = file1.Get(histogram_object_name).GetXaxis().GetXmax()
	Nbins = file1.Get(histogram_object_name).GetSize()
	x_range = np.linspace(xmin, xmax, Nbins+1)
	hist = np.zeros((4, Nbins))
	hist[0,:] = np.linspace(1,Nbins,Nbins)
	hist[1,:] = x_range[0:-1]
	hist[2,:] = x_range[1:]
	hist[3,:] = np.array(file1.Get(histogram_object_name))
	return hist

def read_qkinz_single(isotope, reaction):
	# ** THIS VERSION: Read only the relevant reaction entry (p, d, t) **

	# Function which reads data from Qkinz. 
	# The argument "isotope" should be a string specifying the target, e.g. "184W".
	# The argument "reaction" should be one of the following strings: "p", "d", "t" (for proton, deuteron or triton)
	# Qkinz is available at http://github.com/oslocyclotronlab
	# The Qkinz output should be organized as a set of files named as follows:
	# <isotope>_stripX.txt
	# e.g. 184W_strip0.txt, 184W_strip1.txt, ..., 184W_strip7.txt.
	# The data is returned as nested lists in the following format:
	# [strip1, strip2, ..., strip7],
	# where
	# strip1 = [level1, level2, ..., levelN]
	# where again 
	# level1 = [Excitation energy (keV),     Energy dE-detector (keV),     dEnergy dE-detector (keV),     Energy E-detector (keV),    dEnergy dE-detector (keV),     Total particle energy (keV)]


	list = [] # Allocate nested list to include all data
	for i in range(8):
		filename = "qkinz/%s_strip%d.txt" %(isotope, i)
		infile = open(filename, 'r')
		lines = infile.readlines()
		list_currentstrip = []
		j = 0
		if reaction == "p":
			reactionnumber = 0
		elif reaction == "d":
			reactionnumber = 1
		elif reaction == "t":
			reactionnumber = 2
		else:
			print "read_qkinz: Error in reaction number"
			sys.exit(1)
		readlines_on = False
		firstline = 0
		reactionnumber_counter = 0
		while j < len(lines):
			words = lines[j].split()
			try:
				if str(words[0]) == "Number_of_entries:":
					if reactionnumber == reactionnumber_counter:
						N = int(words[1])
						# print N
						firstline = j + 2
						readlines_on = True
					else:
						reactionnumber_counter += 1 # Not the right reaction yet
			except (ValueError, IndexError): 
				pass
			if readlines_on and j == firstline:
				print "read_qkinz_single consistency check: I am reading a table with", N, "lines."
				for j in range(firstline, firstline+N):
					# Read content of line:
					# Excitation energy (keV):     Energy dE-detector (keV):     dEnergy dE-detector (keV):     Energy E-detector (keV):    dEnergy dE-detector (keV):     Total particle energy (keV): 
					words = lines[j].split()
					list_currentline = [
						float(words[0]),
						float(words[1]),
						float(words[2]),
						float(words[3]),
						float(words[4]),
						float(words[5])						
					]
					list_currentstrip.append(list_currentline)
				break
				# END loop over lines

				

			j += 1
		# End loop over current strip
		list.append(list_currentstrip)
		infile.close()
	# End loop over strips
	return list


# print(read_qkinz_single("184W", "p"))
# sys.exit(0)

def calibrate_2points(peaks_file1, peaks_file2, target1, target2, reaction1, reaction2, excitation_number1, excitation_number2, g0_front = 2.5, g0_back = 5):
	# This function takes two peaks and two matching Qkinz calculated levels, and outputs
	# gains and shifts. 

	# Read measured peak values from files
	infile1 = open(peaks_file1, 'r')
	infile2 = open(peaks_file2, 'r')
	lines1 = infile1.readlines()
	lines2 = infile2.readlines()

	E1_m_list = []
	dE1_m_list = []
	E2_m_list = []
	dE2_m_list = []
	for i in range(64):
		words1 = lines1[2*i+2].split()
		E1_m_list.append(float(words1[1]))
		dE1_m_list.append(float(words1[2]))
		words2 = lines2[2*i+2].split()
		E2_m_list.append(float(words2[1]))
		dE2_m_list.append(float(words2[2]))


	# Get the relevant physical values from Qkinz:
	target1_all_levels = read_qkinz_single(target1, reaction1)
	target2_all_levels = read_qkinz_single(target2, reaction2)

	# Calculate gain and shift for each detector
	# (We use the first measurement to set the gain, and the second to set the shift)
	g_front_list = []
	g_back_list = []
	s_front_list = []
	s_back_list = []
	for i in range(64):
		# All front strips are traversed for each back strip, so e.g. i=2 means front=2, back=0, i=8 means front=0, back=1, etc.

		# Get relevant Qkinz values:
		dE1_c = target1_all_levels[i%8][excitation_number1][1]
		dE2_c = target2_all_levels[i%8][excitation_number2][1]
		E1_c = target1_all_levels[i%8][excitation_number1][3]
		E2_c = target2_all_levels[i%8][excitation_number2][3]

		# Get measurements for current detector
		E1_m = E1_m_list[i]
		dE1_m = dE1_m_list[i]
		E2_m = E2_m_list[i]
		dE2_m = dE2_m_list[i]

		# Calculate front and back gain
		g_front = g0_front*(dE2_c - dE1_c)/(dE2_m - dE1_m)
		g_back = g0_back*(E2_c - E1_c)/(E2_m - E1_m)
		s_front = dE2_c - g_front*dE2_m/g0_front
		s_back = E2_c - g_back*E2_m/g0_back

		# Save results
		g_front_list.append(g_front)
		g_back_list.append(g_back) 
		s_front_list.append(s_front)
		s_back_list.append(s_back) 



	# Print results in a nice way for pasting into gainshifts file
	print "E gain back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_back_list[8*i + 0], g_back_list[8*i + 1], g_back_list[8*i + 2], g_back_list[8*i + 3], g_back_list[8*i + 4], g_back_list[8*i + 5], g_back_list[8*i + 6], g_back_list[8*i + 7])
	print "E gain front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_front_list[8*i + 0], g_front_list[8*i + 1], g_front_list[8*i + 2], g_front_list[8*i + 3], g_front_list[8*i + 4], g_front_list[8*i + 5], g_front_list[8*i + 6], g_front_list[8*i + 7])
	print "E shift back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_back_list[8*i + 0], s_back_list[8*i + 1], s_back_list[8*i + 2], s_back_list[8*i + 3], s_back_list[8*i + 4], s_back_list[8*i + 5], s_back_list[8*i + 6], s_back_list[8*i + 7])
	print "E shift front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_front_list[8*i + 0], s_front_list[8*i + 1], s_front_list[8*i + 2], s_front_list[8*i + 3], s_front_list[8*i + 4], s_front_list[8*i + 5], s_front_list[8*i + 6], s_front_list[8*i + 7])


def line_up_2points(peaks_file1, peaks_file2, peak1_coordinates, peak2_coordinates, g0_front = 1, g0_back = 1):
	# This function takes two peaks and two coordinates of peaks to calibrate to. This function is only useful for lining up detector spectra to peaks 
	# - to calibrate to physical values, use calibrate_2points() instead.

	# Read measured peak values from files
	infile1 = open(peaks_file1, 'r')
	infile2 = open(peaks_file2, 'r')
	lines1 = infile1.readlines()
	lines2 = infile2.readlines()

	E1_m_list = []
	dE1_m_list = []
	E2_m_list = []
	dE2_m_list = []
	for i in range(64):
		words1 = lines1[2*i+2].split()
		E1_m_list.append(float(words1[0]))
		dE1_m_list.append(float(words1[1]))
		words2 = lines2[2*i+2].split()
		E2_m_list.append(float(words2[0]))
		dE2_m_list.append(float(words2[1]))
		# print E1_m_list[-1], dE1_m_list[-1], E2_m_list[-1], dE2_m_list[-1]


	# Set desired peak values:
	E1_c = peak1_coordinates[0]
	dE1_c = peak1_coordinates[1]
	E2_c = peak2_coordinates[0]
	dE2_c = peak2_coordinates[1]

	# Calculate gain and shift for each detector
	# (We use the first measurement to set the gain, and the second to set the shift)
	g_front_list = []
	g_back_list = []
	s_front_list = []
	s_back_list = []
	for i in range(64):
		# All front strips are traversed for each back strip, so e.g. i=2 means front=2, back=0, i=8 means front=0, back=1, etc.

		# Get measurements for current detector
		E1_m = E1_m_list[i]
		dE1_m = dE1_m_list[i]
		E2_m = E2_m_list[i]
		dE2_m = dE2_m_list[i]

		# Calculate front and back gain
		g_front = g0_front*(dE2_c - dE1_c)/(dE2_m - dE1_m)
		g_back = g0_back*(E2_c - E1_c)/(E2_m - E1_m)
		s_front = dE2_c - g_front*dE2_m/g0_front
		s_back = E2_c - g_back*E2_m/g0_back

		# Save results
		g_front_list.append(g_front)
		g_back_list.append(g_back) 
		s_front_list.append(s_front)
		s_back_list.append(s_back) 



	# Print results in a nice way for pasting into gainshifts file
	print "E gain back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_back_list[8*i + 0], g_back_list[8*i + 1], g_back_list[8*i + 2], g_back_list[8*i + 3], g_back_list[8*i + 4], g_back_list[8*i + 5], g_back_list[8*i + 6], g_back_list[8*i + 7])
	print "E gain front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_front_list[8*i + 0], g_front_list[8*i + 1], g_front_list[8*i + 2], g_front_list[8*i + 3], g_front_list[8*i + 4], g_front_list[8*i + 5], g_front_list[8*i + 6], g_front_list[8*i + 7])
	print "E shift back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_back_list[8*i + 0], s_back_list[8*i + 1], s_back_list[8*i + 2], s_back_list[8*i + 3], s_back_list[8*i + 4], s_back_list[8*i + 5], s_back_list[8*i + 6], s_back_list[8*i + 7])
	print "E shift front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_front_list[8*i + 0], s_front_list[8*i + 1], s_front_list[8*i + 2], s_front_list[8*i + 3], s_front_list[8*i + 4], s_front_list[8*i + 5], s_front_list[8*i + 6], s_front_list[8*i + 7])


def calibrate_2points_from_list_summed_backs(peaks1, peaks2, target1, target2, reaction1, reaction2, excitation_number1, excitation_number2, g0_front = 1, g0_back = 1):
	# This function takes two peaks and two matching Qkinz calculated levels, and outputs
	# gains and shifts. 
	# For this version the peaks should be given as lists on the form [[x_coordinates],[y_coordinates]], with only eight entries since the same is used for all different back strips

	E1_m_list = peaks1[0]
	dE1_m_list = peaks1[1]
	E2_m_list = peaks2[0]
	dE2_m_list = peaks2[1]


	# Get the relevant physical values from Qkinz:
	target1_all_levels = read_qkinz_single(target1, reaction1)
	target2_all_levels = read_qkinz_single(target2, reaction2)

	# Calculate gain and shift for each detector
	# (We use the first measurement to set the gain, and the second to set the shift)
	g_front_list = []
	g_back_list = []
	s_front_list = []
	s_back_list = []
	for i in range(64):
		# All front strips are traversed for each back strip, so e.g. i=2 means front=2, back=0, i=8 means front=0, back=1, etc.

		# Get relevant Qkinz values:
		dE1_c = target1_all_levels[i%8][excitation_number1][1]
		dE2_c = target2_all_levels[i%8][excitation_number2][1]
		E1_c = target1_all_levels[i%8][excitation_number1][3]
		E2_c = target2_all_levels[i%8][excitation_number2][3]

		# Get measurements for current detector
		E1_m = E1_m_list[i%8]
		dE1_m = dE1_m_list[i%8]
		E2_m = E2_m_list[i%8]
		dE2_m = dE2_m_list[i%8]

		# Calculate front and back gain
		g_front = g0_front*(dE2_c - dE1_c)/(dE2_m - dE1_m)
		g_back = g0_back*(E2_c - E1_c)/(E2_m - E1_m)
		s_front = dE2_c - g_front*dE2_m/g0_front
		s_back = E2_c - g_back*E2_m/g0_back

		# Save results
		g_front_list.append(g_front)
		g_back_list.append(g_back) 
		s_front_list.append(s_front)
		s_back_list.append(s_back) 

	# Print results in a nice way for pasting into gainshifts file
	print "E gain back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_back_list[8*i + 0], g_back_list[8*i + 1], g_back_list[8*i + 2], g_back_list[8*i + 3], g_back_list[8*i + 4], g_back_list[8*i + 5], g_back_list[8*i + 6], g_back_list[8*i + 7])
	print "E gain front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_front_list[8*i + 0], g_front_list[8*i + 1], g_front_list[8*i + 2], g_front_list[8*i + 3], g_front_list[8*i + 4], g_front_list[8*i + 5], g_front_list[8*i + 6], g_front_list[8*i + 7])
	print "E shift back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_back_list[8*i + 0], s_back_list[8*i + 1], s_back_list[8*i + 2], s_back_list[8*i + 3], s_back_list[8*i + 4], s_back_list[8*i + 5], s_back_list[8*i + 6], s_back_list[8*i + 7])
	print "E shift front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_front_list[8*i + 0], s_front_list[8*i + 1], s_front_list[8*i + 2], s_front_list[8*i + 3], s_front_list[8*i + 4], s_front_list[8*i + 5], s_front_list[8*i + 6], s_front_list[8*i + 7])


def calibrate_2points_simple(E1_measured, E2_measured, E1_reference, E2_reference, g0 = 1):
	# This function takes two measured points and two reference points and calculates gain and shift to apply to the measured points to get them to hit the reference points.

	g = g0*(E2_reference - E1_reference)/(E2_measured - E1_measured)
	s = E2_reference - g*E2_measured/g0

	return g, s

def evaluate_point(E, g, s):
	# Takes an energy value and a gain+shift, returns gained and shifted energy value
	return s + E*g


def calibrate_2points_from_list_hybrid(peaks_file1, peaks2, target1, target2, reaction1, reaction2, excitation_number1, excitation_number2, g0_front = 1, g0_back = 1):
	# This function takes two peaks and two matching Qkinz calculated levels, and outputs
	# gains and shifts. 
	# For this version, the first peak should be a .csv file from peakfinder, and the other should be given as lists on the form [[x_coordinates],[y_coordinates]], with only eight entries since the same is used for all different back strips

	E2_m_list = peaks2[0]
	dE2_m_list = peaks2[1]

	# Read measured peak values from files
	infile1 = open(peaks_file1, 'r')
	lines1 = infile1.readlines()

	E1_m_list = []
	dE1_m_list = []
	for i in range(64):
		words1 = lines1[2*i+2].split()
		E1_m_list.append(float(words1[0]))
		dE1_m_list.append(float(words1[1]))

	# Get the relevant physical values from Qkinz:
	target1_all_levels = read_qkinz_single(target1, reaction1)
	target2_all_levels = read_qkinz_single(target2, reaction2)

	# Calculate gain and shift for each detector
	# (We use the first measurement to set the gain, and the second to set the shift)
	g_front_list = []
	g_back_list = []
	s_front_list = []
	s_back_list = []
	for i in range(64):
		# All front strips are traversed for each back strip, so e.g. i=2 means front=2, back=0, i=8 means front=0, back=1, etc.

		# Get relevant Qkinz values:
		dE1_c = target1_all_levels[i%8][excitation_number1][1]
		dE2_c = target2_all_levels[i%8][excitation_number2][1]
		E1_c = target1_all_levels[i%8][excitation_number1][3]
		E2_c = target2_all_levels[i%8][excitation_number2][3]

		# Get measurements for current detector
		E1_m = E1_m_list[i]
		dE1_m = dE1_m_list[i]
		E2_m = E2_m_list[i%8]
		dE2_m = dE2_m_list[i%8]

		# Calculate front and back gain
		g_front = g0_front*(dE2_c - dE1_c)/(dE2_m - dE1_m)
		g_back = g0_back*(E2_c - E1_c)/(E2_m - E1_m)
		s_front = dE2_c - g_front*dE2_m/g0_front
		s_back = E2_c - g_back*E2_m/g0_back

		# Save results
		g_front_list.append(g_front)
		g_back_list.append(g_back) 
		s_front_list.append(s_front)
		s_back_list.append(s_back) 

	# Print results in a nice way for pasting into gainshifts file
	print "E gain back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_back_list[8*i + 0], g_back_list[8*i + 1], g_back_list[8*i + 2], g_back_list[8*i + 3], g_back_list[8*i + 4], g_back_list[8*i + 5], g_back_list[8*i + 6], g_back_list[8*i + 7])
	print "E gain front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(g_front_list[8*i + 0], g_front_list[8*i + 1], g_front_list[8*i + 2], g_front_list[8*i + 3], g_front_list[8*i + 4], g_front_list[8*i + 5], g_front_list[8*i + 6], g_front_list[8*i + 7])
	print "E shift back:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_back_list[8*i + 0], s_back_list[8*i + 1], s_back_list[8*i + 2], s_back_list[8*i + 3], s_back_list[8*i + 4], s_back_list[8*i + 5], s_back_list[8*i + 6], s_back_list[8*i + 7])
	print "E shift front:"
	for i in range(8):
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(s_front_list[8*i + 0], s_front_list[8*i + 1], s_front_list[8*i + 2], s_front_list[8*i + 3], s_front_list[8*i + 4], s_front_list[8*i + 5], s_front_list[8*i + 6], s_front_list[8*i + 7])


def make_composite_gainshifts_ede(original_gainshift_file, relative_gainshift_file, total_gainshift_file):
	# This function takes a gainshifts file and another file with relative gains and shifts.
	# These should be linked in the following way:
	# The original gainshifts file (e.g. plain-gain) has been used to sort the data,
	# and this sorting has been used to recalibrate the spectra (e.g. lining up detectors), 
	# giving rise to a new set of gains and shifts which are relative to the original ones.
	# This function puts these together and writes a new, composite gainshifts file.

	original_gainshifts_ = open(original_gainshift_file, 'r') # Should be a proper gainshifts.dat file
	original_gainshifts = original_gainshifts_.readlines() 
	relative_gainshifts_ = open(relative_gainshift_file, 'r') # Can be a proper gainshifts file, but can also contain only the first blocks giving gains and shifts for E and dE. Note that the number of lines should be proper, i.e. the lines representing NaI and Ge gains must be present, but they can be empty
	relative_gainshifts = relative_gainshifts_.readlines() 
	total_gainshifts = open(total_gainshift_file, 'w') # File to write results to. Will be formatted as a proper gainshifts.dat


	# Read data
	original_g_E = np.zeros((8,8))
	for i in range(8):
		iline = i
		currentline = original_gainshifts[iline].split()
		original_g_E[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	original_g_dE = np.zeros((8,8))
	for i in range(8):
		iline = i+10
		currentline = original_gainshifts[iline].split()
		original_g_dE[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	original_s_E = np.zeros((8,8))
	for i in range(8):
		iline = i+28
		currentline = original_gainshifts[iline].split()
		original_s_E[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	original_s_dE = np.zeros((8,8))
	for i in range(8):
		iline = i+37
		currentline = original_gainshifts[iline].split()
		original_s_dE[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	relative_g_E = np.zeros((8,8))
	for i in range(8):
		iline = i
		currentline = relative_gainshifts[iline].split()
		relative_g_E[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	relative_g_dE = np.zeros((8,8))
	for i in range(8):
		iline = i+10
		currentline = relative_gainshifts[iline].split()
		relative_g_dE[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	relative_s_E = np.zeros((8,8))
	for i in range(8):
		iline = i+28
		currentline = relative_gainshifts[iline].split()
		relative_s_E[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))
	relative_s_dE = np.zeros((8,8))
	for i in range(8):
		iline = i+37
		currentline = relative_gainshifts[iline].split()
		relative_s_dE[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]), float(currentline[7]))

	# Calculate cumulative gains and shifts
	total_g_E = relative_g_E*original_g_E
	total_s_E = relative_g_E*original_s_E + relative_s_E
	total_g_dE = relative_g_dE*original_g_dE
	total_s_dE = relative_g_dE*original_s_dE + relative_s_dE

	# Write results to file
	# E and dE gains
	for i in range(8):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_g_E[i,0], total_g_E[i,1], total_g_E[i,2], total_g_E[i,3], total_g_E[i,4], total_g_E[i,5], total_g_E[i,6], total_g_E[i,7]))
	total_gainshifts.write("\n")
	total_gainshifts.write("\n")
	for i in range(8):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_g_dE[i,0], total_g_dE[i,1], total_g_dE[i,2], total_g_dE[i,3], total_g_dE[i,4], total_g_dE[i,5], total_g_dE[i,6], total_g_dE[i,7]))
	total_gainshifts.write("\n")
	total_gainshifts.write("\n")
	# Copy Germanium and NaI lines from original
	for iline in range(20,27):
		# print original_gainshifts[iline]
		total_gainshifts.write(original_gainshifts[iline])
	total_gainshifts.write("\n")

	# E and dE shifts
	for i in range(8):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_s_E[i,0], total_s_E[i,1], total_s_E[i,2], total_s_E[i,3], total_s_E[i,4], total_s_E[i,5], total_s_E[i,6], total_s_E[i,7]))
	total_gainshifts.write("\n")
	for i in range(8):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_s_dE[i,0], total_s_dE[i,1], total_s_dE[i,2], total_s_dE[i,3], total_s_dE[i,4], total_s_dE[i,5], total_s_dE[i,6], total_s_dE[i,7]))
	total_gainshifts.write("\n")

	# Copy last part
	for iline in range(46,len(original_gainshifts)):
		total_gainshifts.write(original_gainshifts[iline])

	# Close files
	original_gainshifts_.close()
	relative_gainshifts_.close()
	total_gainshifts.close()

def make_composite_gainshifts_nai(original_gainshift_file, relative_gainshift_file, total_gainshift_file):
	# This function takes a gainshifts file and another file with relative gains and shifts.
	# These should be linked in the following way:
	# The original gainshifts file (e.g. plain-gain) has been used to sort the data,
	# and this sorting has been used to recalibrate the spectra (e.g. lining up detectors), 
	# giving rise to a new set of gains and shifts which are relative to the original ones.
	# This function puts these together and writes a new, composite gainshifts file.

	original_gainshifts_ = open(original_gainshift_file, 'r') # Should be a proper gainshifts.dat file
	original_gainshifts = original_gainshifts_.readlines() 
	relative_gainshifts_ = open(relative_gainshift_file, 'r') # Can be a proper gainshifts file, but can also contain only the first blocks giving gains and shifts for NaI. Note that the number of lines should be proper.
	relative_gainshifts = relative_gainshifts_.readlines() 
	total_gainshifts = open(total_gainshift_file, 'w') # File to write results to. Will be formatted as a proper gainshifts.dat


	# Read data
	original_g_NaI = np.zeros((5,7))
	for i in range(4):
		iline = i+22
		currentline = original_gainshifts[iline].split()
		original_g_NaI[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]))
	currentline = original_gainshifts[26].split()
	original_g_NaI[4,0:4] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]))

	original_s_NaI = np.zeros((5,7))
	for i in range(4):
		iline = i+48
		currentline = original_gainshifts[iline].split()
		original_s_NaI[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]))
	currentline = original_gainshifts[52].split()
	original_s_NaI[4,0:4] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]))

	relative_g_NaI = np.zeros((5,7))
	for i in range(4):
		iline = i+22
		currentline = relative_gainshifts[iline].split()
		relative_g_NaI[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]))
	currentline = relative_gainshifts[26].split()
	relative_g_NaI[4,0:4] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]))

	relative_s_NaI = np.zeros((5,7))
	for i in range(4):
		iline = i+48
		currentline = relative_gainshifts[iline].split()
		relative_s_NaI[i,:] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]), 
							 float(currentline[4]), float(currentline[5]), float(currentline[6]))
	currentline = relative_gainshifts[52].split()
	relative_s_NaI[4,0:4] = (float(currentline[0]), float(currentline[1]), float(currentline[2]), float(currentline[3]))

	# Calculate cumulative gains and shifts
	total_g_NaI = relative_g_NaI*original_g_NaI
	total_s_NaI = relative_g_NaI*original_s_NaI + relative_s_NaI
	
	# == Write results to file ==
	# Copy E and dE lines from original
	for iline in range(22):
		# print original_gainshifts[iline]
		total_gainshifts.write(original_gainshifts[iline])

	# Write NaI gains
	for i in range(4):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_g_NaI[i,0], total_g_NaI[i,1], total_g_NaI[i,2], total_g_NaI[i,3], total_g_NaI[i,4], total_g_NaI[i,5], total_g_NaI[i,6]))
	total_gainshifts.write("%.5f %.5f %.5f %.5f\n" %(total_g_NaI[4,0], total_g_NaI[4,1], total_g_NaI[4,2], total_g_NaI[4,3]))

	# Copy E and dE lines from original
	for iline in range(27,47):
		# print original_gainshifts[iline]
		total_gainshifts.write(original_gainshifts[iline])
	total_gainshifts.write("\n")

	# NaI shifts
	for i in range(4):
		total_gainshifts.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %(total_s_NaI[i,0], total_s_NaI[i,1], total_s_NaI[i,2], total_s_NaI[i,3], total_s_NaI[i,4], total_s_NaI[i,5], total_s_NaI[i,6]))
	total_gainshifts.write("%.5f %.5f %.5f %.5f\n" %(total_s_NaI[4,0], total_s_NaI[4,1], total_s_NaI[4,2], total_s_NaI[4,3]))

	# Copy last part
	for iline in range(53,len(original_gainshifts)):
		total_gainshifts.write(original_gainshifts[iline])

	# Close files
	original_gainshifts_.close()
	relative_gainshifts_.close()
	total_gainshifts.close()



def read_ede_histograms(filename):
	infile = open(filename, 'r')
	lines = infile.readlines()
	
	i_lines = 0
	
	back = []
	i_back = 0
	# Read histograms for back detectors (assuming 64 different spectra, i.e. gated on front strips)
	while i_back < 64:
		Nbins = int(lines[i_lines].split()[-1])
		back_current = np.zeros((4,Nbins))
		i_lines += 2
		for i in range(Nbins):
			words = lines[i+i_lines].split()
			back_current[:,i] = (float(words[0]), float(words[1]), float(words[2]), float(words[3]))
			# print back_current[:,i]
		back.append(back_current)
		i_lines += Nbins
		i_back += 1
	
	# Repeat for front detectors
	front = []
	i_front = 0
	while i_front < 64:
		Nbins = int(lines[i_lines].split()[-1])
		front_current = np.zeros((4,Nbins))
		i_lines += 2
		for i in range(Nbins):
			words = lines[i+i_lines].split()
			front_current[:,i] = (float(words[0]), float(words[1]), float(words[2]), float(words[3]))
			# print front_current[:,i]
		front.append(front_current)
		i_lines += Nbins
		i_front += 1
	
	return front, back

def read_nai_histograms(filename):
	infile = open(filename, 'r')
	lines = infile.readlines()
	
	
	i_lines = 0
	
	spectra = []
	i_spectra = 0
	# Read histograms for NaI detectors
	while i_spectra < 32:
		Nbins = int(lines[i_lines].split()[-1])
		spectra_current = np.zeros((4,Nbins))
		i_lines += 2
		for i in range(Nbins):
			words = lines[i+i_lines].split()
			spectra_current[:,i] = (float(words[0]), float(words[1]), float(words[2]), float(words[3]))
			# print back_current[:,i]
		spectra.append(spectra_current)
		i_lines += Nbins
		i_spectra += 1

	return spectra

def rebin_data(data, factor):
	# Factors of 2, please (or identity)
	if factor % 2 != 0 and factor != 1:
		print "Function rebin only handles powers of two (and identity)"
		sys.exit(1)
	if factor == 1:
		return data
	else:
		Nbins = len(data[0,:])
		data_new = np.zeros((4, int(Nbins/factor)))
		for i in range(int(Nbins/factor)):
			data_new[:,i] = (i, data[1,factor*i], data[2,factor*(i+1)-1], np.sum(data[3,factor*i:factor*(i+1)]))
		return data_new


def line_up_ede(spectra, i_current,  i_reference=0, Nbins=[1000], binshift=[0], rebin=1, error_option=1, shift_width=50, show_plots=True, normalize_spectra=True):
	# Function to line up spectra using a squared-difference minimization routine.
	#
	# Has the option to use only parts (one or several) of the spectra, and also to rebin (only powers of two for the rebinning factor, please).
	# There are three different options for how the error function to minimize is calculated, using either the difference,
	# the difference between the derivatives or the difference between logarithms of values.

	# The middle of each bin is used as the energy value

	# Get number of calibration regions
	Nregions = len(Nbins)

	if rebin != 1:
		spectrum_ref = rebin_data(spectra[i_reference], rebin)
		spectrum_cur = rebin_data(spectra[i_current], rebin)
	else:
		spectrum_ref = spectra[i_reference]
		spectrum_cur = spectra[i_current]
	# Set reference
	reference = np.array([(spectrum_ref[1,:]+spectrum_ref[2,:])/2,spectrum_ref[3,:]])
	# Choose which detector to line up
	current = np.array([(spectrum_cur[1,:]+spectrum_cur[2,:])/2,spectrum_cur[3,:]])

	# Should probabaly normalize the area of current to reference
	if normalize_spectra:
		A_ref = np.sum(reference[1,binshift[0]:Nbins[0]+binshift[0]])
		A_cur = np.sum(current[1,binshift[0]:Nbins[0]+binshift[0]])
		if Nregions > 1:
			for i in range(1,Nregions):
				A_ref += np.sum(reference[1,binshift[i]:Nbins[i]+binshift[i]])
				A_cur += np.sum(current[1,binshift[i]:Nbins[i]+binshift[i]])
		current[1,:] *= A_ref/A_cur
	
	# Idea: Use an array of indices which is (0,1,2,3..) for the reference and slightly stretched and shifted for the calibration spectrum.
	reference_indices = np.linspace(binshift[0],Nbins[0]-1+binshift[0], Nbins[0])
	if Nregions > 1:
		for i in range(1,Nregions):
			reference_indices = np.append(reference_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	# print Nregions
	# print reference_indices
	
	def error_function(coeff):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		current_indices = np.linspace(binshift[0], Nbins[0]-1+binshift[0], Nbins[0]) # These will be corrected by the minimization
		if Nregions > 1:
			for i in range(1,Nregions):
				current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
		# current_indices = coeff[2]*current_indices**2 + coeff[1]*current_indices + coeff[0] # Use this instead if you want second-order calibration
		current_indices = coeff[1]*current_indices + coeff[0]
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				error = np.sum(np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2))
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[1,reference_indices.astype(int)])-np.diff(current[1,current_indices.astype(int)]),2))
			# Option 3: The error is the sum of squared differences between logarithms
			elif error_option == 3:
				minval=1e-20 # Clip array to avoid taking log of zero
				error = np.sum(np.power(np.log(reference[1,reference_indices.astype(int)].clip(min=minval))-np.log(current[1,current_indices.astype(int)].clip(min=minval)),2))
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)

		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*coeff[1]**2
		return error
	
	# Run minimization. Scan over starting values for shift parameter to ensure we find something optimal.
	min_fval = 1e10
	coefficients = [-1,-1]
	for c0 in range(-shift_width,shift_width,2):
		x0 = [c0, 1.0]
		result = minimize(error_function, x0, method='Nelder-Mead')
		if result.success and result.fun <= min_fval:
			min_fval = result.fun
			# c0_best = c0
			coefficients = result.x
	# print c0_best, coefficients
	
	
	current_indices = np.linspace(binshift[0], Nbins[0]-1+binshift[0], Nbins[0]) # These will be corrected by the minimization
	if Nregions > 1:
		for i in range(1,Nregions):
			current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	# current_indices = coefficients[2]*current_indices**2 + coefficients[1]*current_indices + coefficients[0]
	current_indices = coefficients[1]*current_indices + coefficients[0]
	current_indices[current_indices<0] = 0 # Remove negative values, important since Python interprets negative indices as starting from the top
	
	
	if show_plots:
		# Plot to check that it looks reasonable
		plt.hold('on')
	
		plt.figure(1)
		plt.step(reference[0,reference_indices.astype(int)], reference[1,reference_indices.astype(int)])
		plt.step(current[0,reference_indices.astype(int)], current[1,current_indices.astype(int)])
		plt.step(current[0,reference_indices.astype(int)], current[1,reference_indices.astype(int)])
		plt.title('Spectra')
		plt.legend(['Reference', 'Current, aligned', 'Current, not aligned'], fontsize=12, loc='best')
		if error_option == 3:
			plt.yscale('log')
	
		plt.figure(2)
		plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2), color='g')
		plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,reference_indices.astype(int)],2), color='r')
		plt.title('Squared error')
		plt.legend(['reference minus aligned', 'reference minus not aligned'], fontsize=12, loc='best')
		if error_option == 3:
			plt.yscale('log')
		plt.show()

	# print coefficients
	# print reference_indices, current_indices
	# print reference_indices[-1]*coefficients[1] - current_indices[-1]
	# binsize = reference[0,1]-reference[0,0]
	# print "Binsize =", binsize
	# print "binsize * shift =", binsize*coefficients[0]
	# print reference[0,reference_indices[-1]], reference[0,current_indices[-1].astype(int)]
	# print reference[0,reference_indices.astype(int)][-1]*coefficients[1]
	# print reference[0,current_indices.astype(int)]
	# print reference[0,reference_indices[-1].astype(int)]*coefficients[1]-reference[0,current_indices[-1].astype(int)]

	# Calculate gain and shift coefficients the stupid but (hopefully) failsafe way:
	g = (reference[0,reference_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)])/(reference[0,current_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)])
	s = reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)]-g*reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)]
	# print g,s

	# print reference[0,reference_indices.astype(int)]-(reference[0,current_indices.astype(int)]*g+s)


	# plt.figure(3)
	# plt.plot(reference[0,:], reference[1,:])
	# plt.plot((current[0,:]*g+s), current[1,:])
	# plt.plot(current[0,:], current[1,:])
	# plt.show()

	return g, s



def line_up_ede_better_plotting(spectra, i_current,  i_reference=0, Nbins=[1000], binshift=[0], rebin=1, error_option=1, shift_width=50, show_plots=True, normalize_spectra=True):
	# Function to line up spectra using a squared-difference minimization routine.
	#
	# Has the option to use only parts (one or several) of the spectra, and also to rebin (only powers of two for the rebinning factor, please).
	# There are three different options for how the error function to minimize is calculated, using either the difference,
	# the difference between the derivatives or the difference between logarithms of values.

	# The middle of each bin is used as the energy value

	# Get number of calibration regions
	Nregions = len(Nbins)

	if rebin != 1:
		spectrum_ref = rebin_data(spectra[i_reference], rebin)
		spectrum_cur = rebin_data(spectra[i_current], rebin)
	else:
		spectrum_ref = spectra[i_reference]
		spectrum_cur = spectra[i_current]
	# Set reference
	reference = np.array([(spectrum_ref[1,:]+spectrum_ref[2,:])/2,spectrum_ref[3,:]])
	# Choose which detector to line up
	current = np.array([(spectrum_cur[1,:]+spectrum_cur[2,:])/2,spectrum_cur[3,:]])

	# Should probabaly normalize the area of current to reference
	if normalize_spectra:
		A_ref = np.sum(reference[1,binshift[0]:Nbins[0]+binshift[0]])
		A_cur = np.sum(current[1,binshift[0]:Nbins[0]+binshift[0]])
		if Nregions > 1:
			for i in range(1,Nregions):
				A_ref += np.sum(reference[1,binshift[i]:Nbins[i]+binshift[i]])
				A_cur += np.sum(current[1,binshift[i]:Nbins[i]+binshift[i]])
		current[1,:] *= A_ref/A_cur
	
	# Idea: Use an array of indices which is (0,1,2,3..) for the reference and slightly stretched and shifted for the calibration spectrum.
	reference_indices = np.linspace(binshift[0],Nbins[0]-1+binshift[0], Nbins[0])
	if Nregions > 1:
		for i in range(1,Nregions):
			reference_indices = np.append(reference_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	# print Nregions
	# print reference_indices
	
	def error_function(coeff):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		current_indices = np.linspace(binshift[0], Nbins[0]-1+binshift[0], Nbins[0]) # These will be corrected by the minimization
		if Nregions > 1:
			for i in range(1,Nregions):
				current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
		# current_indices = coeff[2]*current_indices**2 + coeff[1]*current_indices + coeff[0] # Use this instead if you want second-order calibration
		current_indices = coeff[1]*current_indices + coeff[0]
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				error = np.sum(np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2))
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[1,reference_indices.astype(int)])-np.diff(current[1,current_indices.astype(int)]),2))
			# Option 3: The error is the sum of squared differences between logarithms
			elif error_option == 3:
				minval=1e-20 # Clip array to avoid taking log of zero
				error = np.sum(np.power(np.log(reference[1,reference_indices.astype(int)].clip(min=minval))-np.log(current[1,current_indices.astype(int)].clip(min=minval)),2))
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)

		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*coeff[1]**2
		return error
	
	# Run minimization. Scan over starting values for shift parameter to ensure we find something optimal.
	min_fval = 1e10
	coefficients = [-1,-1]
	for c0 in range(-shift_width,shift_width,2):
		x0 = [c0, 1.0]
		result = minimize(error_function, x0, method='Nelder-Mead')
		if result.success and result.fun <= min_fval:
			min_fval = result.fun
			# c0_best = c0
			coefficients = result.x
	# print c0_best, coefficients
	
	
	current_indices = np.linspace(binshift[0], Nbins[0]-1+binshift[0], Nbins[0]) # These will be corrected by the minimization
	if Nregions > 1:
		for i in range(1,Nregions):
			current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	# current_indices = coefficients[2]*current_indices**2 + coefficients[1]*current_indices + coefficients[0]
	current_indices = coefficients[1]*current_indices + coefficients[0]
	current_indices[current_indices<0] = 0 # Remove negative values, important since Python interprets negative indices as starting from the top
	
	
	# if show_plots:
	# 	# Plot to check that it looks reasonable
	# 	plt.hold('on')
	
	# 	plt.figure(1)
	# 	plt.step(reference[0,reference_indices.astype(int)], reference[1,reference_indices.astype(int)])
	# 	plt.step(current[0,reference_indices.astype(int)], current[1,current_indices.astype(int)])
	# 	plt.step(current[0,reference_indices.astype(int)], current[1,reference_indices.astype(int)])
	# 	plt.title('Spectra')
	# 	plt.legend(['Reference', 'Current, aligned', 'Current, not aligned'], fontsize=12, loc='best')
	# 	if error_option == 3:
	# 		plt.yscale('log')
	
	# 	plt.figure(2)
	# 	plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2), color='g')
	# 	plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,reference_indices.astype(int)],2), color='r')
	# 	plt.title('Squared error')
	# 	plt.legend(['reference minus aligned', 'reference minus not aligned'], fontsize=12, loc='best')
	# 	if error_option == 3:
	# 		plt.yscale('log')
	# 	plt.show()




	# print coefficients
	# print reference_indices, current_indices
	# print reference_indices[-1]*coefficients[1] - current_indices[-1]
	# binsize = reference[0,1]-reference[0,0]
	# print "Binsize =", binsize
	# print "binsize * shift =", binsize*coefficients[0]
	# print reference[0,reference_indices[-1]], reference[0,current_indices[-1].astype(int)]
	# print reference[0,reference_indices.astype(int)][-1]*coefficients[1]
	# print reference[0,current_indices.astype(int)]
	# print reference[0,reference_indices[-1].astype(int)]*coefficients[1]-reference[0,current_indices[-1].astype(int)]

	# Calculate gain and shift coefficients the stupid but (hopefully) failsafe way:
	g = (reference[0,reference_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)])/(reference[0,current_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)])
	s = reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)]-g*reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)]
	# print g,s

	# print reference[0,reference_indices.astype(int)]-(reference[0,current_indices.astype(int)]*g+s)

	if show_plots:
		spectra = [rebin_data(spectrum, rebin) for spectrum in spectra]
		N = len(spectra[0][0,:])
		binsize = spectra[i_reference][2,0]-spectra[i_reference][1,0]
		plt.figure(1)
		plt.hold('on')
		plt.step((spectra[i_reference][1,0:N]+spectra[i_reference][2,0:N])/2 , spectra[i_reference][3,0:N], label='spectrum %d, reference'%i_reference, color='r')
		plt.step(((spectra[i_current][1,0:N]+spectra[i_current][2,0:N])/2)*g+s , spectra[i_current][3,0:N], label='spectrum %d, aligned'%i_current, color='b')
		plt.step(((spectra[i_current][1,0:N]+spectra[i_current][2,0:N])/2) , spectra[i_current][3,0:N], label='spectrum %d, not aligned'%i_current, color='c')
		plt.xlim([spectra[i_reference][1,0], spectra[i_reference][1,max((binshift[0]+Nbins[0])*1.2, N*0.95)]])
		if error_option == 3:
			plt.yscale('log')
		plt.title('Line-up calibration check')
		for i in range(Nregions):
			plt.axvline(x=binshift[i]*binsize, color='black', linewidth=2)
			plt.axvline(x=(binshift[i]+Nbins[i])*binsize, color='black', linewidth=2)
			plt.text(x=(binshift[i]+Nbins[i]/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region %d" %i, horizontalalignment='center')
		plt.legend(loc='best')
		plt.show()


	# plt.figure(3)
	# plt.plot(reference[0,:], reference[1,:])
	# plt.plot((current[0,:]*g+s), current[1,:])
	# plt.plot(current[0,:], current[1,:])
	# plt.show()

	return g, s





def line_up_nai(spectra, i_current,  i_reference=0, Nbins=[1000], binshift=[0], rebin=1, error_option=1, shift_width=50, show_plots=True, show_error_plots=False, normalize_spectra=True):
	# Function to line up spectra using a squared-difference minimization routine.
	#
	# Has the option to use only parts of the spectra, and also to rebin in case of low statistics.
	# There are three different options for how the error function to minimize is calculated, using either the difference,
	# the difference between the derivatives or a combination. The last one might need some tweaking to get the relative weight 
	# correct.

	# We use the middle of each bin as the energy value

	# Get number of calibration regions
	Nregions = len(Nbins)

	if rebin != 1:
		spectrum_ref = rebin_data(spectra[i_reference], rebin)
		spectrum_cur = rebin_data(spectra[i_current], rebin)
	else:
		spectrum_ref = spectra[i_reference]
		spectrum_cur = spectra[i_current]
	# Set reference
	reference = np.array([(spectrum_ref[1,:]+spectrum_ref[2,:])/2,spectrum_ref[3,:]])
	# Nbins = 1000 # Number of bins to care about
	# Choose which detector to line up
	current = np.array([(spectrum_cur[1,:]+spectrum_cur[2,:])/2,spectrum_cur[3,:]])

	# Should probabaly normalize the area of current to reference
	if normalize_spectra:
		A_ref = np.sum(reference[1,binshift[0]:Nbins[0]+binshift[0]])
		A_cur = np.sum(current[1,binshift[0]:Nbins[0]+binshift[0]])
		if Nregions > 1:
			for i in range(1,Nregions):
				A_ref += np.sum(reference[1,binshift[i]:Nbins[i]+binshift[i]])
				A_cur += np.sum(current[1,binshift[i]:Nbins[i]+binshift[i]])
		current[1,:] *= A_ref/A_cur
	
	# Idea: Use an array of indices which is (0,1,2,3..) for the reference and slightly stretched and shifted for the calibration spectrum.
	reference_indices = np.linspace(binshift[0],Nbins[0]-1+binshift[0], Nbins[0])
	if Nregions > 1:
		for i in range(1,Nregions):
			reference_indices = np.append(reference_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	
	def error_function(coeff):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		current_indices = np.linspace(binshift[0], Nbins[0]-1+binshift[0], Nbins[0]) # These will be corrected by the minimization
		if Nregions > 1:
			for i in range(1,Nregions):
				current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
		# current_indices = coeff[2]*current_indices**2 + coeff[1]*current_indices + coeff[0] # Use this instead if you want second-order calibration
		current_indices = coeff[1]*current_indices + coeff[0]
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				error = np.sum(np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2))
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[1,reference_indices.astype(int)])-np.diff(current[1,current_indices.astype(int)]),2))
			# Option 3: The error is taken to be the sum of squared differences between logs. 
			elif error_option == 3:
				minval = 1e-20 # Clip array to avoid taking log of zero
				error = np.sum(np.power(np.log(reference[1,reference_indices.astype(int)].clip(min=minval))-np.log(current[1,current_indices.astype(int)].clip(min=minval)),2))
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)

		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*coeff[1]**2
		return error
	
	# Run minimization. Scan over starting values for shift parameter to ensure we find something optimal.
	min_fval = 1e10
	coefficients = [-1,-1]
	for c0 in range(-shift_width,shift_width,2):
		x0 = [c0, 1.0]
		result = minimize(error_function, x0, method='Nelder-Mead')
		if result.success and result.fun <= min_fval:
			min_fval = result.fun
			# c0_best = c0
			coefficients = result.x
	# print c0_best, coefficients
	
	
	current_indices = np.linspace(binshift[0],Nbins[0]-1+binshift[0], Nbins[0]) # This will be corrected by the minimization
	if Nregions > 1:
		for i in range(1,Nregions):
			current_indices = np.append(current_indices, np.linspace(binshift[i],Nbins[i]-1+binshift[i], Nbins[i]))
	# current_indices = coefficients[2]*current_indices**2 + coefficients[1]*current_indices + coefficients[0]
	current_indices = coefficients[1]*current_indices + coefficients[0]
	current_indices[current_indices<0] = 0 # Remove negative values, important since Python interprets negative indices as starting from the top
	
	


	# print coefficients
	# print reference_indices, current_indices
	# print reference_indices[-1]*coefficients[1] - current_indices[-1]
	# binsize = reference[0,1]-reference[0,0]
	# print "Binsize =", binsize
	# print "binsize * shift =", binsize*coefficients[0]
	# print reference[0,reference_indices[-1]], reference[0,current_indices[-1].astype(int)]
	# print reference[0,reference_indices.astype(int)][-1]*coefficients[1]
	# print reference[0,current_indices.astype(int)]
	# print reference[0,reference_indices[-1].astype(int)]*coefficients[1]-reference[0,current_indices[-1].astype(int)]

	# Calculate gain and shift coefficients the stupid but (hopefully) failsafe way:
	g = (reference[0,reference_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)])/(reference[0,current_indices[int(np.sum(Nbins)*3/4)].astype(int)]-reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)])
	s = reference[0,reference_indices[int(np.sum(Nbins)*1/4)].astype(int)]-g*reference[0,current_indices[int(np.sum(Nbins)*1/4)].astype(int)]
	# print g,s

	# print reference[0,reference_indices.astype(int)]-(reference[0,current_indices.astype(int)]*g+s)


	# plt.figure(3)
	# plt.plot(reference[0,:], reference[1,:])
	# plt.plot((current[0,:]*g+s), current[1,:])
	# plt.plot(current[0,:], current[1,:])
	# plt.show()


	# == Plotting ==
	binsize = reference[0,1]-reference[0,0]
	N = len(reference[0,:])
	
	if show_plots:
		# Plot to check that it looks reasonable
		plt.hold('on')
	

		plt.figure(1)
		
		plt.step((spectrum_ref[1,0:N]+spectrum_ref[2,0:N])/2 , spectrum_ref[3,0:N], label='spectrum %d'%i_reference, color="blue")
		# plt.xlim([spectrum_ref[1,0], spectrum_ref[1,max((binshift+Nbins)*1.2, N*0.95)]])
		
		plt.step(((spectrum_cur[1,0:N]+spectrum_cur[2,0:N])/2) , spectrum_cur[3,0:N], label='spectrum %d, not corrected'%i_current, color="darkseagreen")
		plt.step(((spectrum_cur[1,0:N]+spectrum_cur[2,0:N])/2)*g+s , spectrum_cur[3,0:N], label='spectrum %d, corrected'%i_current, color="crimson")
	
		if error_option == 3:
			plt.yscale('log')

		# plt.step(reference[0,reference_indices.astype(int)], reference[1,reference_indices.astype(int)])
		# plt.step(current[0,reference_indices.astype(int)], current[1,current_indices.astype(int)])
		# plt.step(current[0,reference_indices.astype(int)], current[1,reference_indices.astype(int)])
		plt.title('Spectra')
		plt.legend(fontsize=12, loc='best')
		for i in range(Nregions):
			plt.axvline(x=binshift[i]*binsize, color='black', linewidth=2)
			plt.axvline(x=(binshift[i]+Nbins[i])*binsize, color='black', linewidth=2)
			plt.text(x=(binshift[i]+Nbins[i]/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region %d" %i, horizontalalignment='center')
	if show_error_plots:
		plt.figure(2)
		plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,reference_indices.astype(int)],2), color='darkseagreen', label='reference minus uncalibrated')
		plt.step(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2), color='crimson', label='reference minus calibrated')
		plt.title('Squared error')
		plt.legend(fontsize=12, loc='best')

		if error_option == 3:
			plt.yscale('log')

		for i in range(Nregions):
			plt.axvline(x=binshift[i]*binsize, color='black', linewidth=2)
			plt.axvline(x=(binshift[i]+Nbins[i])*binsize, color='black', linewidth=2)
			plt.text(x=(binshift[i]+Nbins[i]/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region %d" %i, horizontalalignment='center')
	if show_plots or show_error_plots:
		plt.show()




	return g, s


def line_up_ede_2D(current, reference, Nbins_x, Nbins_y, Nbins_tot, binshift_x=0, binshift_y=0, error_option=1, shift_width=20, shift_step=4, show_plots=True, normalize_spectra=True, original_Nbins_x=2000, original_Nbins_y=2000, low_value_threshold=0, high_value_threshold=100, axis_limits=[0,40000,0,14000]):
	

	Nbins_tot_x, Nbins_tot_y = Nbins_tot


	# Make linearly spaced axis numbering. This is used later to calculate gains and shifts.
	axis_x = np.linspace(axis_limits[0], axis_limits[1], Nbins_tot_x)
	axis_y = np.linspace(axis_limits[2], axis_limits[3], Nbins_tot_y)
	
	# Make arrays containing indices. One set is kept as reference, the other is shifted and scaled to fit.
	ref_indices_x = np.linspace(binshift_x, Nbins_x-1+binshift_x, Nbins_x)
	ref_indices_y = np.linspace(binshift_y, Nbins_y-1+binshift_y, Nbins_y)
	
	# Normalize area of current to reference
	if normalize_spectra:
		A_ref = np.sum(reference[binshift_x:Nbins_x+binshift_x,binshift_y:Nbins_y+binshift_y])
		A_cur = np.sum(reference[binshift_x:Nbins_x+binshift_x,binshift_y:Nbins_y+binshift_y])
		current[:,:] *= A_ref/A_cur
	
	
	# test = np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))], 2)
	# print test
	
	
	colorbar_min = 0
	colorbar_max = np.log((np.max(reference) + np.max(current))**2)
	
	if show_plots:
		fig, axes = plt.subplots(2,2)
	
		logsquarediffspectra = np.log(np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))], 2))
		logsquaresumspectra = np.log(np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]+current[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))], 2))
	
		
		# plt.subplot(2,2,1)
		axes.flat[0].set_title('$\log((\mathrm{current} + \mathrm{reference})^2)$, not aligned')
		axes.flat[0].imshow(logsquaresumspectra, origin='lower', vmin=colorbar_min, vmax=colorbar_max)
		# plt.colorbar()
	
		# plt.subplot(2,2,3)
		axes.flat[2].set_title('$\log((\mathrm{current} - \mathrm{reference})^2)$, not aligned, inside value threshold')
		# Set all elements outside threshold limits to zero
		logsquarediffspectra[np.logical_or(logsquarediffspectra<low_value_threshold, logsquarediffspectra>high_value_threshold)] = np.nan
		axes.flat[2].imshow(logsquarediffspectra, origin='lower', vmin=colorbar_min, vmax=colorbar_max)
		# extent=[axis_x[binshift_x], axis_x[binshift_x+Nbins_x],
		#	  axis_y[binshift_y], axis_y[binshift_y+Nbins_y]])
	
		# plt.colorbar()
	
	
	def error_function(coeff):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		
		# The coeff[] list is assumed as [shift_x, gain_x, shift_y, gain_y].
		cur_indices_x = coeff[1]*ref_indices_x + coeff[0]
		cur_indices_y = coeff[3]*ref_indices_y + coeff[2]
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				squarediffspectra = np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))], 2)
				squarediffspectra = squarediffspectra[np.logical_and(squarediffspectra > np.exp(low_value_threshold), squarediffspectra < np.exp(high_value_threshold))] # cut away low values
				error = np.sum(np.log(squarediffspectra))
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[ref_indices_x.astype(int),ref_indices_y.astype(int)])-np.diff(current[ref_indices_x.astype(int),ref_indices_y.astype(int)]),2))
			elif error_option == 3:
				error = np.sum(np.power(np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))])-np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))]),2))
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)
	
		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*(coeff[1]**2 + coeff[3]**2)
		return error
	
	# min_fval = 1e10
	# coefficients = [-1,-1,-1,-1]
	# for c0_x in range(-shift_width,shift_width,shift_step):
	# 	for c0_y in range(-shift_width,shift_width,shift_step):
	# 		x0 = [c0_x, 1.0, c0_y, 1.0]
	# 		result = minimize(error_function, x0, method='Nelder-Mead')
	# 		if result.success and result.fun <= min_fval:
	# 			min_fval = result.fun
	# 			# c0_best = c0
	# 			coefficients = result.x
	# print coefficients
	
	def error_function_y(coeff_y, coeff_x0, coeff_x1):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		
		# The coeff[] list is assumed as [shift_x, gain_x, shift_y, gain_y].
		cur_indices_x = coeff_x1*ref_indices_x + coeff_x0
		cur_indices_y = coeff_y[1]*ref_indices_y + coeff_y[0]
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				squarediffspectra = np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))], 2)
				squarediffspectra = squarediffspectra[np.logical_and(squarediffspectra > np.exp(low_value_threshold), squarediffspectra < np.exp(high_value_threshold))] # cut away low values
				error = np.sum(np.log(squarediffspectra))
				# error = np.sum(squarediffspectra)
				# error = np.sum(np.power(np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))])-np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.	astype(int))]),2)
					# *np.power(reference_indices,1.5) # Multiply by the index, to give some kind of increasing weight to data further to the right
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[ref_indices_x.astype(int),ref_indices_y.astype(int)])-np.diff(current[ref_indices_x.astype(int),ref_indices_y.astype(int)]),2))
			elif error_option == 3:
				squarediffspectra = np.power(np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))])-np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))]), 2)
				squarediffspectra = squarediffspectra[np.logical_and(squarediffspectra > np.exp(low_value_threshold), squarediffspectra < np.exp(high_value_threshold))] # cut away low values
				error = np.sum(squarediffspectra)
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)
	
		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*(coeff_x1**2 + coeff_y[1]**2)
		return error
	
	def error_function_x(coeff_x, coeff_y0, coeff_y1):
		# The squared difference of the spectra, given calibration coefficients. This will be minimized.
		
		# The coeff[] list is assumed as [shift_x, gain_x, shift_y, gain_y].
		cur_indices_x = coeff_x[1]*ref_indices_x + coeff_x[0]
		cur_indices_y = coeff_y1*ref_indices_y + coeff_y0
		# print coeff
		try:
			# Option 1: The error is taken to be the sum of squared differences. 
			if error_option == 1:
				squarediffspectra = np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))], 2)
				squarediffspectra = squarediffspectra[np.logical_and(squarediffspectra > np.exp(low_value_threshold), squarediffspectra < np.exp(high_value_threshold))] # cut away low values
				error = np.sum(np.log(squarediffspectra))
				# error = np.sum(np.power(np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))])-np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))]),2)
						# *np.power(reference_indices,1.5) # Multiply by the index, to give some kind of increasing weight to data further to the right
					
			# Option 2: The error is the sum of squared differences between derivatives			
			elif error_option == 2:
				error = np.sum(np.power(np.diff(reference[ref_indices_x.astype(int),ref_indices_y.astype(int)])-np.diff(current[ref_indices_x.astype(int),ref_indices_y.astype(int)]),2))
			elif error_option == 3:
				squarediffspectra = np.power(np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))])-np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))]), 2)
				squarediffspectra = squarediffspectra[np.logical_and(squarediffspectra > np.exp(low_value_threshold), squarediffspectra < np.exp(high_value_threshold))] # cut away low values
				error = np.sum(squarediffspectra)
			else:
				print "Error in function error_function: Bad value for error_option."
				sys.exit(1)
	
		except IndexError: # Sometimes the minimizer chooses coefficients which overstep the bounds of the array
			error = 1e10*(coeff_x[1]**2 + coeff_y1**2)
		return error
	
	min_fval = 1e10
	# First minimize in x direction
	coefficients_x = [-1,-1]
	for c0_x in range(-shift_width,shift_width,shift_step):
		x0 = [c0_x, 1.0]
		y0 = [0, 1.0]
		result = minimize(error_function_x, x0, args=(y0[0], y0[1]), method='Nelder-Mead')
		if result.success and result.fun <= min_fval:
			min_fval = result.fun
			# c0_best = c0
			coefficients_x = result.x
	# min_fval_x = min_fval
	# print coefficients_x
	# Then minimize in y direction
	coefficients_y = [-1, -1]
	for c0_y in range(-shift_width,shift_width,shift_step):
		x0 = coefficients_x
		y0 = [c0_y, 1.0]
		result = minimize(error_function_y, y0, args=(x0[0], x0[1]), method='Nelder-Mead')
		if result.success and result.fun <= min_fval:
			min_fval = result.fun
			# c0_best = c0
			coefficients_y = result.x
	# min_fval_y = min_fval
	# print coefficients_y
	# Rerun for x
	result = minimize(error_function_x, coefficients_x, args=(coefficients_y[0], coefficients_y[1]), method='Nelder-Mead')
	if result.success and result.fun < min_fval:
		print "Updated min_fval in x-direction afterburner (from, to):", min_fval, result.fun
		min_fval = result.fun
		# c0_best = c0
		print "Corresponding coefficients (from, to):", coefficients_x, result.x
		coefficients_x = result.x
	# Rerun for y
	result = minimize(error_function_y, coefficients_y, args=(coefficients_x[0], coefficients_x[1]), method='Nelder-Mead')
	if result.success and result.fun < min_fval:
		print "Updated min_fval in y-direction afterburner (from, to):", min_fval, result.fun
		min_fval = result.fun
		# c0_best = c0
		print "Corresponding coefficients (from, to):", coefficients_y, result.x
		coefficients_y = result.x
	# Finally try all four coefficients at once
	coefficients = [coefficients_x[0], coefficients_x[1], coefficients_y[0], coefficients_y[1]]
	result = minimize(error_function, coefficients, method='Nelder-Mead')
	if result.success and result.fun < min_fval:
		min_fval = result.fun
		# c0_best = c0
		coefficients = result.x
		print "Updated min_fval in four-way afterburner minimization (from, to):", min_fval, result.fun
		print "Corresponding coefficients (from, to):", coefficients, result.x


	


	# Make array of best shift indices
	cur_indices_x = np.copy(ref_indices_x)
	cur_indices_x = coefficients[1]*cur_indices_x + coefficients[0]
	cur_indices_x[cur_indices_x<0] = 0 # Remove negative values, important since Python interprets negative indices as starting from the top
	cur_indices_y = np.copy(ref_indices_y)
	cur_indices_y = coefficients[3]*cur_indices_y + coefficients[2]
	cur_indices_y[cur_indices_y<0] = 0
	
	
	if show_plots:
		# Plot to check that it looks reasonable
		# plt.hold('on')

		# plt.figure(2)

		logsquarediffspectra = np.log(np.power(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]-current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))], 2))
		logsquaresumspectra = 2*np.log(reference[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]+current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))])

		# plt.subplot(2,2,2)
		axes.flat[1].set_title('$\log((\mathrm{current} + \mathrm{reference})^2)$, aligned')
		axes.flat[1].imshow(logsquaresumspectra, origin='lower', vmin=colorbar_min, vmax=colorbar_max)
		# plt.colorbar()

		# plt.subplot(2,2,4)
		# Set all elements outside threshold limits to zero
		axes.flat[3].set_title('$\log((\mathrm{current} - \mathrm{reference})^2)$, aligned, inside value threshold')
		logsquarediffspectra[np.logical_or(logsquarediffspectra<low_value_threshold, logsquarediffspectra>high_value_threshold)] = np.nan
		im = axes.flat[3].imshow(logsquarediffspectra, origin='lower', vmin=colorbar_min, vmax=colorbar_max)
		# plt.colorbar()

		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		fig.colorbar(im, cax=cbar_ax)

	
		#extent=[axis_x[binshift_x], axis_x[binshift_x+Nbins_x], axis_y[binshift_y], axis_y[binshift_y+Nbins_y]])
		# plt.imshow(np.log(current[np.meshgrid(ref_indices_x.astype(int),ref_indices_y.astype(int))]))
		# plt.imshow(np.log(current[np.meshgrid(cur_indices_x.astype(int),cur_indices_y.astype(int))]))
		# plt.legend(['Reference', 'Current, aligned'], fontsize=12, loc='upper left')
	
		# plt.figure(2)
		# plt.plot(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,current_indices.astype(int)],2), color='g')
		# plt.plot(reference[0,reference_indices.astype(int)], np.power(reference[1,reference_indices.astype(int)]-current[1,reference_indices.astype(int)],2), color='r')
		# plt.title('Squared error')
		# plt.legend(['reference minus calibrated', 'reference minus uncalibrated'], fontsize=12, loc='upper left')
	plt.show()


	# Calculate gains and shifts
	g_x = (axis_x[ref_indices_x[int(Nbins_x*3/4)].astype(int)]-axis_x[ref_indices_x[int(Nbins_x*1/4)].astype(int)])/(axis_x[cur_indices_x[int(Nbins_x*3/4)].astype(int)]-axis_x[cur_indices_x[int(Nbins_x*1/4)].astype(int)])
	s_x = axis_x[ref_indices_x[int(Nbins_x*1/4)].astype(int)]-g_x*axis_x[cur_indices_x[int(Nbins_x*1/4)].astype(int)]
	g_y = (axis_y[ref_indices_y[int(Nbins_y*3/4)].astype(int)]-axis_y[ref_indices_y[int(Nbins_y*1/4)].astype(int)])/(axis_y[cur_indices_y[int(Nbins_y*3/4)].astype(int)]-axis_y[cur_indices_y[int(Nbins_y*1/4)].astype(int)])
	s_y = axis_y[ref_indices_y[int(Nbins_y*1/4)].astype(int)]-g_y*axis_y[cur_indices_y[int(Nbins_y*1/4)].astype(int)]


	return s_x, g_x, s_y, g_y




def calibrate_ede_strip(spectra, strip, Nbins=[1000], rebin=1, binshift=[0], error_option=1, shift_width=50, show_plots=True, normalize_spectra=True):
	# Function to calibrate a whole strip of front/back detectors. Spectra=front/back, strip=0,1,...,7.
	# 
	# Binshift indicates the number of *indices* to skip before reaching the part of the spectra used for calibration (i.e. squared-error fitting). 
	# It can be translated into a shift on the x-axis of the plot by multiplying with the bin size.
	# Nbins similarly is the number of indices to include in the calibration procedure.
	#
	spectra = [rebin_data(spectrum, rebin) for spectrum in spectra] # Rebin here rather than inside the calibration, to enable the function to plot with the same bin size
	binsize = spectra[strip][2,0]-spectra[0][1,0]
	print "Binsize =", binsize
	N = len(spectra[0][0,:])
	Nregions = len(Nbins)
	if show_plots:
		plt.figure(1)
		plt.hold('on')
		plt.step((spectra[strip][1,0:N]+spectra[strip][2,0:N])/2 , spectra[strip][3,0:N], label='spectrum %d'%strip)
		plt.xlim([spectra[strip][1,0], spectra[strip][1,int(max((binshift[0]+Nbins[0])*1.2, N*0.95))]])
		if error_option == 3:
			plt.yscale('log')
		plt.title('Line-up calibration check')
		for i in range(Nregions):
			plt.axvline(x=binshift[i]*binsize, color='black', linewidth=2)
			plt.axvline(x=(binshift[i]+Nbins[i])*binsize, color='black', linewidth=2)
			plt.text(x=(binshift[i]+Nbins[i]/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region %d" %i, horizontalalignment='center')
	print "Detector #	gain 		shift "
	print "%d			1			0" %strip
	gains = [1]
	shifts = [0]
	for i in range(1,8):
		i = 8*i+strip
		g,s = line_up_ede(spectra, i, strip, Nbins=Nbins, binshift=binshift, rebin=1, error_option=error_option, shift_width=shift_width, show_plots=False, normalize_spectra=normalize_spectra)
		gains.append(g)
		shifts.append(s)
		print "%d 			%f 	%f" %(i,g,s)
		if show_plots:
			plt.step(((spectra[i][1,0:N]+spectra[i][2,0:N])/2)*g+s , spectra[i][3,0:N], label='spectrum %d'%i)
	
	if show_plots:
		plt.legend(fontsize=10)
		plt.show()
	return gains, shifts



def calibrate_ede_all_individual_detector_tune(spectra, Nbins, rebin, binshift, error_option, shift_width, show_plots=True, normalize_spectra=False):
	# Function to calibrate a whole strip of front/back detectors. Spectra=front/back, strip=0,1,...,7.
	# 
	# Binshift indicates the number of *indices* to skip before reaching the part of the spectra used for calibration (i.e. squared-error fitting). 
	# It can be translated into a shift on the x-axis of the plot by multiplying with the bin size.
	# Nbins similarly is the number of indices to include in the calibration procedure.
	#

	# Rebin here rather than inside the calibration, to enable the function to plot with the same bin size
	spectra_tmp = []
	for i_back in range(8):
		for i_front in range(8):
			spectra_tmp.append(rebin_data(spectra[8*i_back + i_front], rebin[i_front][i_back]))
			# print "i_back = %d, i_front = %d, spectrum = %d, rebin = %d" %(i_back, i_front, 8*i_back + i_front, rebin[i_front][i_back])
	spectra = spectra_tmp 

	binsize = spectra[0][2,0]-spectra[0][1,0]
	print "Binsize =", binsize
	N = len(spectra[0][0,:])

	gains = []
	shifts = []
	for i_front in range(8):
		if show_plots:
			plt.figure(i_front)
			plt.hold('on')
			plt.step((spectra[i_front][1,0:N]+spectra[i_front][2,0:N])/2 , spectra[i_front][3,0:N], label='spectrum %d'%i_front)
			plt.xlim([spectra[i_front][1,0], spectra[i_front][1,max((binshift[i_front][1][0]+Nbins[i_front][1][0])*1.2, N*0.95)]])
			if error_option[i_front][1] == 3:
				plt.yscale('log' )
			plt.title('Line-up calibration check, detector ring %d' %i_front)
			Nregions = len(Nbins[i_front][1]) # Assuming all detectors in same strip have same number of regions. Also using only Nbins/binshift from detector 1 for plotting, so might be inaccurate
			for i in range(Nregions):
				plt.axvline(x=binshift[i_front][1][i]*binsize, color='black', linewidth=2)
				plt.axvline(x=(binshift[i_front][1][i]+Nbins[i_front][1][i])*binsize, color='black', linewidth=2)
				plt.text(x=(binshift[i_front][1][i]+Nbins[i_front][1][i]/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region %d" %i, horizontalalignment='center')
		print "Detector #	gain 		shift "
		print "%d			1			0" %i_front
		gains.append(1)
		shifts.append(0)
		for i_back in range(1,8):
			# i_back = 8*i_back+strip
			g,s = line_up_ede_better_plotting(spectra, (8*i_back+i_front), i_front, Nbins=Nbins[i_front][i_back], binshift=binshift[i_front][i_back], rebin=1, error_option=error_option[i_front][i_back], shift_width=shift_width[i_front][i_back], show_plots=False, normalize_spectra=normalize_spectra)
			gains.append(g)
			shifts.append(s)
			print "%d 			%f 	%f" %(8*i_back+i_front,g,s)
			if show_plots:
				plt.step(((spectra[8*i_back+i_front][1,0:N]+spectra[8*i_back+i_front][2,0:N])/2)*g+s, spectra[8*i_back+i_front][3,0:N], label='spectrum %d'%(8*i_back+i_front))
		
		if show_plots:
			plt.legend(fontsize=10)
			plt.show()

	# Print gains and shifts nice for pasting into gainshifts file:
	print "Gains:"
	for i in range(8):
		# print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(gains[8*i + 0], gains[8*i + 1], gains[8*i + 2], gains[8*i + 3], gains[8*i + 4], gains[8*i + 5], gains[8*i + 6], gains[8*i + 7])
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(gains[i + 0*8], gains[i + 1*8], gains[i + 2*8], gains[i + 3*8], gains[i + 4*8], gains[i + 5*8], gains[i + 6*8], gains[i + 7*8])
	print "Shifts:"
	for i in range(8):
		# print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(shifts[8*i + 0], shifts[8*i + 1], shifts[8*i + 2], shifts[8*i + 3], shifts[8*i + 4], shifts[8*i + 5], shifts[8*i + 6], shifts[8*i + 7])
		print "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(shifts[i + 0*8], shifts[i + 1*8], shifts[i + 2*8], shifts[i + 3*8], shifts[i + 4*8], shifts[i + 5*8], shifts[i + 6*8], shifts[i + 7*8])

	return gains, shifts



def calibrate_nai(spectra, indices, Nbins, rebin, binshift, error_option=[1], shift_width=[50], show_plots=True, normalize_spectra=True, print_on=True, logarithmic=True):
	# Function to calibrate all NaI detectors. 
	# Uses spectra[0] as reference.
	#
	# Binshift indicates the number of *indices* to skip before reaching the part of the spectra used for calibration (i.e. squared-error fitting). 
	# It can be translated into a shift on the x-axis of the plot by multiplying with the bin size.
	# Nbins similarly is the number of indices to include in the calibration procedure.
	#
	# spectra = [rebin_data(spectrum, rebin) for spectrum in spectra] # Rebin here rather than inside the calibration, to enable the function to plot with the same bin size
	binsize = spectra[0][2,0]-spectra[0][1,0]
	print "Binsize =", binsize
	N = len(spectra[0][0,:])
	if show_plots:
		plt.figure(1)
		plt.hold('on')
		plt.step((spectra[0][1,0:N]+spectra[0][2,0:N])/2 , spectra[0][3,0:N], label='spectrum 0')
		# plt.xlim([spectra[0][1,0], spectra[0][1,max((binshift[0]+Nbins[0])*1.2, N*0.95)]])
		plt.ylim([0,3*max(spectra[0][3,:])])
		plt.title('Line-up calibration check')
		# plt.axvline(x=binshift*binsize, color='black', linewidth=2)
		# plt.axvline(x=(binshift+Nbins)*binsize, color='black', linewidth=2)
		# plt.text(x=(binshift+Nbins/2)*binsize, y=max(spectra[0][3,:]), s="Calibration region", horizontalalignment='center')
		if logarithmic:
			plt.yscale('log')
	if print_on:
		print "Detector #	gain 		shift "
		print "0			1			0" 
	gains = [1]
	shifts = [0]
	for i in range(1,len(indices)):
		if indices[i]:
			g,s = line_up_nai(spectra, i, Nbins=Nbins[i], binshift=binshift[i], rebin=rebin[i], error_option=error_option[i], shift_width=shift_width[i], show_plots=False, normalize_spectra=normalize_spectra)
			gains.append(g)
			shifts.append(s)
			if print_on:
				print "%d 			%f 	%f" %(i,g,s)
			if show_plots:
				plt.step(((spectra[i][1,0:N]+spectra[i][2,0:N])/2)*g+s , spectra[i][3,0:N], label='spectrum %d'%i)
		else:
			gains.append(0)
			shifts.append(0)
	if show_plots:
		plt.legend(fontsize=8)
		plt.show()
	return gains, shifts


def fit_gauss(spectrum, fit_range, strip_number=-1, show_plots=True):

    from scipy.optimize import curve_fit
    def gauss(x, mu, sigma, C):
        return C*np.exp(-np.power(x-mu,2)/(2*np.power(sigma,2)))
    def log_gauss(x, mu, sigma, C):
        return np.log(C) -np.power(x-mu,2)/(2*np.power(sigma,2))


    # Do boring stuff
    spectrum_x = (spectrum[1,:]+spectrum[2,:])/2
    spectrum_y = spectrum[3,:]
    #print np.where(spectrum_x>fit_range[0])[0][0]
    index_range = [np.where(spectrum_x>fit_range[0])[0][0], np.where(spectrum_x>fit_range[1])[0][0]]
    #print index_range
    # Estimate starting values for parameters
    mu0 = np.mean(fit_range)
    sigma0 = np.sqrt(np.sum(np.power(spectrum_x[index_range[0]:index_range[1]]-mu0,2)))
    C0 = spectrum_y[int(np.mean(index_range))]
    p0 = [mu0, sigma0, C0]

    #print "p0 =", p0

    popt1, pcov1 = curve_fit(gauss, spectrum_x[index_range[0]:index_range[1]], spectrum_y[index_range[0]:index_range[1]], p0)
    popt2, pcov2 = curve_fit(gauss, spectrum_x[index_range[0]:index_range[1]], np.log(spectrum_y[index_range[0]:index_range[1]]), p0)

    #print "popt, normal =", popt1
    #print "popt, log =", popt2

    #plt.figure()
    if show_plots:
    	plt.figure()
    	plt.hold('on')
    	plt.step(spectrum_x, spectrum_y, label="spectrum")
    	if strip_number >= 0:
    		plt.title('Strip number %d' %strip_number)
    	plt.plot(spectrum_x, gauss(spectrum_x, popt1[0], popt1[1], popt1[2]), label="regular fit")
    	plt.plot(spectrum_x, np.exp(gauss(spectrum_x, popt2[0], popt2[1], popt2[2])), label="log fit")
    	plt.yscale('log')
    	plt.legend(loc='best')
    	plt.ylim([1e0, 1e6])
    	plt.show()
    #plt.hold('off')
    
    return popt1, popt2



def print_gainshifts_block_nai(g, s):
	# Simply takes a gain and a shift parameter and prints blocks formatted to fit gainshifts files. 
	# Typically used when spectra have already been aligned, and only physical calibration remains.
	print "Gains:"
	for i in range(4):
		i *= 7
		print "%f %f %f %f %f %f %f" %(g, g, g, g, g, g, g) 		
	print "%f %f %f %f" %(g, g, g, g)
	print "Shifts:"
	for i in range(4):																
		i *= 7
		print "%f %f %f %f %f %f %f" %(s, s, s, s, s, s, s) 		
	print "%f %f %f %f" %(s, s, s, s)
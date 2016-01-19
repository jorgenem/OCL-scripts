from __future__ import division
import numpy as np
import sys

def read_qkinz(isotope, protons_on = True, deuterons_on = False, tritons_on = False):
	# Function which reads data from Qkinz. 
	# The Qkinz output should be organized as a set of files named as follows:
	# <isotope>_stripX.txt
	# e.g. 184W_strip0.txt, 184W_strip1.txt, ..., 184W_strip7.txt.
	# The data is returned as nested lists in the following format:
	# [strip1, strip2, ..., strip7],
	# where
	# strip1 = [protons, deuterons, tritons] (or not all these, depending on settings)
	# where again
	# protons = [level1, level2, ..., levelN]
	# where again 
	# level1 = [Excitation energy (keV),     Energy dE-detector (keV),     dEnergy dE-detector (keV),     Energy E-detector (keV),    dEnergy dE-detector (keV),     Total particle energy (keV)]


	list = [] # Allocate nested list to include all data
	for i in range(8):
		filename = "qkinz/%s_strip%d.txt" %(isotope, i)
		infile = open(filename, 'r')
		lines = infile.readlines()
		list_currentstrip = []
		j = 0
		protons_done = False
		deuterons_done = False
		tritons_done = False
		readlines_on = False
		firstline = 0
		while j < len(lines):
			words = lines[j].split()
			try:
				if str(words[0]) == "Number_of_entries:":
					N = int(words[1])
					print N
					firstline = j + 2
					readlines_on = True
			except (ValueError, IndexError): 
				pass
			if readlines_on and j == firstline:
				list_currentreaction = []
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
					list_currentreaction.append(list_currentline)
				# END loop over lines

				list_currentstrip.append(list_currentreaction)
				j = firstline + N - 1 # Set j to the right value to continue in the outside while loop

				# We just finished reading some species. Check which species it was, and whether we should continue with the next or break because we are done.
				# NB! This currently does not work perfectly for all reaction choices. (Problems if protons are off, I think). Use with caution.
				if protons_on:
					if not protons_done:
						protons_done = True
						if not deuterons_on and not tritons_on:
							break # All species done
					else:
						if deuterons_on:
							if not deuterons_done:
								deuterons_done = True
								if not tritons_on:
									break # All species done
							else:
								if tritons_on:
									if not tritons_done:
										tritons_done = True
									else:
										break # All species done
						else: # Protons are on, but deuterons are off
							if tritons_on:
								if not tritons_done:
									tritons_done = True
								else:
									break # All species done
							else: 
								break # All species done
				else: # Protons are not on
					if deuterons_on:
						if not deuterons_done:
							deuterons_done = True
							if not tritons_on:
								break # All species done
						else:
							if tritons_on:
								if not tritons_done:
									tritons_done = True
									break # All species done
								else:
									break # All species done
					else: # Deuterons are not on
						if tritons_on:
							if not tritons_done:
								tritons_done = True # Superfluous unless another species is added later
								break # All species done
							else:
								break # All species done
						else:
							print "No particles on?"
							break				

				readlines_on = False

			j += 1
		# End loop over current strip
		list.append(list_currentstrip)
		infile.close()
	# End loop over strips
	return list


def read_qkinz_single(isotope, reaction):
	# ** THIS VERSION: Read only the relevant reaction entry (p, d, t) **

	# Function which reads data from Qkinz. 
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
				print "read_qkinz_single DEBUG: I am reading a table with", N, "lines."
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
		print "%f %f %f %f %f %f %f %f" %(g_back_list[8*i + 0], g_back_list[8*i + 1], g_back_list[8*i + 2], g_back_list[8*i + 3], g_back_list[8*i + 4], g_back_list[8*i + 5], g_back_list[8*i + 6], g_back_list[8*i + 7])
	print "E gain front:"
	for i in range(8):
		print "%f %f %f %f %f %f %f %f" %(g_front_list[8*i + 0], g_front_list[8*i + 1], g_front_list[8*i + 2], g_front_list[8*i + 3], g_front_list[8*i + 4], g_front_list[8*i + 5], g_front_list[8*i + 6], g_front_list[8*i + 7])
	print "E shift back:"
	for i in range(8):
		print "%f %f %f %f %f %f %f %f" %(s_back_list[8*i + 0], s_back_list[8*i + 1], s_back_list[8*i + 2], s_back_list[8*i + 3], s_back_list[8*i + 4], s_back_list[8*i + 5], s_back_list[8*i + 6], s_back_list[8*i + 7])
	print "E shift front:"
	for i in range(8):
		print "%f %f %f %f %f %f %f %f" %(s_front_list[8*i + 0], s_front_list[8*i + 1], s_front_list[8*i + 2], s_front_list[8*i + 3], s_front_list[8*i + 4], s_front_list[8*i + 5], s_front_list[8*i + 6], s_front_list[8*i + 7])

calibrate_2points(	"triton_ground_state_peak-second_try-stripped.csv", 
					"triton_first_excited_prominent_peak_stripped.csv",
					"184W", "184W", "t", "t", 0, 5, 2.5, 5)
	

	
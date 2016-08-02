from __future__ import division
import numpy as np

isotope = "184W"
# Assuming protons

# Suboptimal formatting of Qkinz table makes this script suboptimal too...

coefficients = []
for i in range(8):
	filename = "%s_strip%d.txt" %(isotope, i)
	infile = open(filename, 'r')
	lines = infile.readlines()
	words = lines[34].split()
	coefficients.append([words[4], words[7], words[10]])

for coefficients_line in coefficients:
	print coefficients_line, "\\"
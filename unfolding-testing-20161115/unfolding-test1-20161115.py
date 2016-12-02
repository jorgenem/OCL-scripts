from __future__ import division

import numpy as np 
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/OCL-scripts')
import pyma_v1 as pyma

import unfolding_test1_f2py as unfolding
print unfolding.__doc__

# Import raw mama matrix
raw, calib, Eg_range, Ex_range = pyma.read_mama('alfna-20160518.m')
R, tmp1, tmp2, tmp3 = pyma.read_mama('response-20161122.m')

raw_reformat = np.zeros((4096,2048))
raw_reformat[0:raw.shape[0],0:raw.shape[1]] = raw
# raw_reformat[0:raw.T.shape[0],0:raw.T.shape[1]] = raw.T # Test with transposition
# raw_reformat[10,100] = 100
print raw_reformat.max()
plt.matshow(raw_reformat,origin='lower')
plt.show()

plt.matshow(R,origin='lower')
plt.show()

# sys.exit(0)

# Run unfolding
# rawmat = np.zeros((4096,2048))
calib = calib[0:2]
print "calib = ",calib
unfmat = unfolding.unfoldit(raw_reformat,calib,R,iversion=2)

print unfmat.max()
plt.matshow(unfmat,origin='lower')
plt.show()
#!/usr/bin/python -i

from membrainrunner import *

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = 20
framecount_force_exact = True
location = ''
execfile('locations.py')

picklename = 'pkl.gr2d.v531.PIP2-PIP2.s4-sim-trestles-md.part0004.xtc.12000-16000-10.pkl'
mset = unpickle(pickles+picklename)

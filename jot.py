#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

mset = unpickle(pickles+'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl')

plt.imshow(std(mset.surf,axis=0).T,interpolation='nearest',origin='lower');plt.show()

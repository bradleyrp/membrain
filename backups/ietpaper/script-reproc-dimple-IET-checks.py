#!/usr/bin/python -i

interact = True
from membrainrunner import *
execfile('locations.py')

structpkldir = 'backup-2014.01.09-enth-review-pickles-exo70-pickles-transpose-error/'

mset = unpickle(pickles+structpkldir+'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl')
dat = mean(mset.surf,axis=0)
plt.imshow(array(mset.surf[50]).T,interpolation='nearest',origin='lower');plt.show()
plt.imshow(array(mset.surf[20]).T,interpolation='nearest',origin='lower');plt.show()


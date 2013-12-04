#!/usr/bin/python -i

from membrainrunner import *
location = 'dark'
execfile('locations.py')
execfile('plotter.py')
mset = unpickle(pickles+'pkl.postproc-dimple-fit.membrane-v614.md.part0002.skip10.pkl')
mset.calculate_undulation_spectrum(removeavg=0,redundant=1)
mset.analyze_undulations(redundant=1)
plotter_undulations_summary(mset,qmagfilter=[0.01,0.5],imagefile=pickles+'undulation-spectrum.png')

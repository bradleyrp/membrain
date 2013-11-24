#!/usr/bin/python -i

from membrainrunner import *
location = 'dark'
execfile('locations.py')

#mset = unpickle('pkl.avgstruct.membrane-v509.md.part0033.skip10.pbcmol.pkl')
mset = unpickle(pickles+'pkl.structures.membrane-v530.md.part0004.skip100.pkl')
mset.calculate_undulation_spectrum(removeavg=1,redundant=1,whichframes=list(where([min([i for j in mset.surf[k] for i in j])>-20. for k in range(len(mset.surf))])[0]))
mset.analyze_undulations(redundant=1)
plotter_undulations_summary(mset,qmagfilter=[.1,.8],zoom=True)

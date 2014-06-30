#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---load the standard header defintions
execfile('header-meso.py')

simcode = 'v2013'
analyses_names = [simcode+'-'+(meso_expt_toc[simcode])['parameter_name']+'-'+str(i) 
	for i in (meso_expt_toc[simcode])['parameter_sweep']]

#---set the interpolation length scale in native simulation units
lenscale = 1.0

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---interpolate structure and induced curvature for mesoscale simulations if not available
lenscale = 0.5
for a in analyses_names:
	for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
	if 'mset' in globals(): del mset
	mset = unpickle(pickles+'pkl.structures.meso.'+(analysis_descriptors[a])['testname']+'.pkl')
	status('status: checking for structure pickle, '+(analysis_descriptors[a])['testname'])
	if mset == None:
		mset = MembraneSet()
		status('status: computing structure, '+(analysis_descriptors[a])['testname'])
		c0sraw = array(mset.load_points_vtu(locate,extra_props='induced_cur',
			start=start,end=end,nbase=nbase,lenscale=lenscale,prefix='EQUIB-conf-'))[:,0]
		mset.surfacer()
		c0s = mset.surfacer_general(c0sraw)
		result_data = MembraneData('c0map')
		result_data.data = array(array(c0s))
		for i in analysis_descriptors[a]: 
			result_data.addnote([i,(analysis_descriptors[a])[i]])
		mset.store.append(result_data)
		pickledump(mset,'pkl.structures.meso.'+(analysis_descriptors[a])['testname']+'.pkl',
			directory=pickles)
	


#!/usr/bin/python

#---PRECOMPUTE MESOSCALE MODEL HEIGHTS AND IMPOSED CURVATURE
#---Pairs well with script-coupling.py with the meso_precomp options

from membrainrunner import *
execfile('locations.py')

import scipy.ndimage
from scipy.signal import hilbert

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
analysis_descriptors = {
	'v2002-t3':
		{'simtype':'meso',
		'shortname':r'meso(iso)',
		'testname':'v2002-t3',
		'locate':
			'/store-delta/compbio/mesoscale-v2002/t3-anis-22-c0-0.05/run1-size-sweep/rep-0/equilibrate/',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':True,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t4':
		{'simtype':'meso',
		'shortname':'meso(bare)',
		'testname':'v2002-t4',
		'locate':\
			'/home/rpb/worker/repo-membrane/mesoscale-v2002/t4-bare-22/run1-size-sweep/rep-0/equilibrate/',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':False,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t2':
		{'simtype':'meso',
		'shortname':r'meso(aniso)',
		'testname':'v2002-t2',
		'locate':
			'/store-delta/compbio/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':True,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t1':
		{'simtype':'meso',
		'shortname':'meso(bare)',
		'testname':'v2002-t1',
		'locate':\
			'/store-delta/compbio/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':False,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v614':
		{'simtype':'md',
		'shortname':r'$4\times$ENTH(MD)',
		'testname':'v614',
		'locate':'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl',
		'start':None,
		'end':None,
		'nbase':None,
		'hascurv':True,
		'hypo':[0.05,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True},
	'v009':
		{'simtype':'meso',
		'shortname':'meso(bare)',
		'testname':'v009',
		'locate':\
			'/home/rpb/compbio-alt/meso-v009-stationary-field-30x30-bare/serial-code/equilibrate/',
		'start':400,
		'end':700,
		'nbase':30,
		'hascurv':False,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True}}

analyses_names = ['v2002-t4','v2002-t3','v2002-t1','v2002-t2','v009'][-1:]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
lenscale = 1.0
for a in analyses_names:
	for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
	if 'mset' in globals(): del mset
	mset = MembraneSet()
	if simtype == 'meso':
		c0sraw = array(mset.load_points_vtu(locate,extra_props='induced_cur',
			start=start,end=end,nbase=nbase,lenscale=lenscale))[:,0]
		mset.surfacer()
		c0s = mset.surfacer_general(c0sraw)
	result_data = MembraneData('c0map')
	result_data.data = array(array(c0s))
	for i in analysis_descriptors[a]: 
		result_data.addnote([i,(analysis_descriptors[a])[i]])
	mset.store.append(result_data)
	pickledump(mset,'pkl.structures.meso.'+(analysis_descriptors[a])['testname']+'.pkl',
		directory=pickles)
	


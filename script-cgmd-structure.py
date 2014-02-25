#!/usr/bin/python

logfile,interact,debugmode = [None,False,None]
from membrainrunner import *
execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = None

#---parameters
rounder = 20.0

#---selections
sel_cgmd_surfacer = ['name PO4 or name POG','name C2A']
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc'},
	'v700-500000-600000-200':
		{'sysname':'membrane-v700','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u1-lonestar-longrun/md.part0009.500000-700000-200.xtc',
		'timeslice':[500000,600000,200]},
	'v701-60000-160000-200':
		{'sysname':'membrane-v701','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc'},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc'},
	'v550-4000000-500000-160':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc'},
	'v550-300000-400000-200':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc'},
		}	
		
analysis_names = ['v701-60000-160000-200','v614-120000-220000-200','v700-500000-600000-200',
	'v612-75000-175000-200','v550-4000000-5000000-160','v550-300000-400000-200'][-1:]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over analysis questions
for aname in analysis_names:
	#---details
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
	#---loop over trajectory files
	for traj in trajfile:
		mset = MembraneSet()
		#---load the trajectory
		print 'status: accessing '+specname_pickle(sysname,traj)
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='cgmd')
		checktime()
		#---average structure calculation
		mset.identify_monolayers(director)
		#---infer the timeslice from the XTC name if not specified
		if 'timeslice' in analysis_descriptors[aname].keys():
			print 'warning: requested timeslice '+str(timeslice)
		else:
			if len(trajsel.split('.')[-2].split('-')) == 1:
				tslicepos = -3
				subset_name = trajsel.split('.')[-2]
			else: tslicepos = -2
			timeslice = [int(i) for i in trajsel.split('.')[tslicepos].split('-')]
		if protein_select == None:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,timeslice=timeslice)
		else:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
				protein_selection=protein_select,timeslice=timeslice)
		mset.calculate_undulations()
		#---save the data
		pickledump(mset,'pkl.structures.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)
		if erase_when_finished:
			del mset
		checktime()


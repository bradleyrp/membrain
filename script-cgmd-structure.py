#!/usr/bin/python

logfile,interact,debugmode = [None,False,None]
from membrainrunner import *
execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---parameters
rounder = 20.0

#---load the standard header defintions
execfile('header-cgmd.py')
		
analysis_names = [
	'v701-60000-160000-200',
	'v614-120000-220000-200',
	'v700-500000-600000-200',
	'v612-75000-175000-200',
	'v550-4000000-5000000-160',
	'v550-300000-400000-200',
	'v614-40000-140000-200',
	'v612-10000-80000-200',
	'v616-210000-310000-200',
	][-1:]

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


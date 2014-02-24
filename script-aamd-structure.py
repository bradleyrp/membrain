#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = None

#---parameters
rounder = 4.0

#---selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector_aamd_symmetric = 'name P'
selector_aamd_asymmetric = '(name P and not resname CHL1)'
selector_aamd_asymmetric = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
residues_aamd_symmetric = ['DOPC','DOPS','PI2P']
residues_aamd_asymmetric = ['DOPC','DOPS','DOPE','POPC','P35P','PI2P']
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-atomP',
		'director':director_aamd_symmetric,'selector':selector_aamd_symmetric,'protein_select':None,
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.atomP.xtc',
		'timeslice':[30000,80000,100]}}
analysis_names = ['v511-30000-80000-100']

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
		basename = traj.split('/')[-1][:-4]
		#---revised basename to include step-part because sometimes the time gets reset
		basename = "-".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])[:-4]
		print 'status: accessing '+basename
		starttime = time.time()
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---average structure calculation
		mset.identify_monolayers(director)
		if protein_select == None:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
				timeslice=timeslice,thick=True)
		else:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
				protein_selection=protein_select,timeslice=timeslice,thick=True)
		mset.calculate_undulations()
		#---save the data
		pickledump(mset,'pkl.structures.'+sysname+'.'+basename[:11]+'.'+str(timeslice[0])+'-'+
			str(timeslice[1])+'-'+str(timeslice[2])+'.pkl',directory=pickles)
		if erase_when_finished:
			del mset
		checktime()


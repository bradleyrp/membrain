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

#---possible analyses
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614',
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'timeslice':[120000,220000,200]},
	'v510-40000-90000-1000':
		{'sysname':'membrane-v510',
		'director':director_symmetric,'selector':selector_aamd,'protein_select':None,
		'trajsel':'u1-lonestar-longrun/md.part0009.500000-700000-200.xtc',
		'timeslice':[40000,90000,1000]}}
analysis_names = ['v510-40000-90000-1000']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time(); print 'start'
#---loop over analysis questions
for aname in analysis_names:
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	#---file lookup
	if type(trajsel) == slice:
		trajfile = trajectories[systems.index(sysname)][trajsel]
	elif type(trajsel) == str:
		pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajsel))
		for fname in trajectories[systems.index(sysname)]:
			if pat.match(fname):
				trajfile = [fname]
	elif type(trajsel) == list:
		trajfile = []
		for trajfilename in trajsel:
			pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajfilename))
			for fname in trajectories[systems.index(sysname)]:
				if pat.match(fname):
					trajfile.append(fname)
	print 'trajectories: '+str(trajfile)
	#---loop over trajectory files
	for traj in trajfile:
		mset = MembraneSet()
		#---Load the trajectory
		gro = structures[systems.index(sysname)]
		basename = traj.split('/')[-1][:-4]
		#---revised basename to include step-part because sometimes the time gets reset
		basename = "-".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])[:-4]
		sel_surfacer = sel_aamd_surfacer
		print 'Accessing '+basename+'.'
		starttime = time.time()
		mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
		print 'time = '+str(1./60*(time.time()-starttime))+' minutes.'
		#---Average structure calculation
		mset.identify_monolayers(director,startframeno=0)
		if protein_select == None:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,timeslice=timeslice)
		else:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
				protein_selection=protein_select,timeslice=timeslice)
		mset.calculate_undulations()
		#---Save the data
		pickledump(mset,'pkl.structures.'+sysname+'.'+basename[:11]+'.'+str(timeslice[0])+'-'+
			str(timeslice[1])+'-'+str(timeslice[2])+'.pkl',directory=pickles)
		if erase_when_finished:
			del mset
		print 'time = '+str(1./60*(time.time()-starttime))+' minutes.'


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

#---analysis plan
analysis_plan = slice(None,None)
analysis_descriptors = [
	('membrane-v614',director_cgmd,selector_cgmd,
		None,'s9-lonestar/md.part0004.120000-220000-200.xtc',[120000,220000,200])]
analyses = [analysis_descriptors[i] for i in [0]]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time(); print 'start'
#---loop over analysis questions
for ad in analyses:
	(sysname,director,selector,protein_select,trajsel,timeslice) = ad
	#---file lookup
	if type(trajsel) == slice:
		trajfile = trajectories[systems.index(sysname)][trajsel]
	elif type(trajsel) == str:
		pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajsel))
		for fname in trajectories[systems.index(sysname)]:
			if pat.match(fname):
				trajfile = [fname]
	elif type(traj) == list:
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
		sel_surfacer = sel_cgmd_surfacer
		print 'Accessing '+basename+'.'
		starttime = time.time()
		mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='cgmd')
		print 'time = '+str(1./60*(time.time()-starttime))+' minutes.'
		#---Average structure calculation
		mset.identify_monolayers(director,startframeno=0)
		if protein_select == None:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,timeslice=timeslice)
		else:
			mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
				protein_selection=protein_select,timeslice=timeslice)
		mset.calculate_undulation_spectrum(removeavg=0,redundant=1)
		mset.analyze_undulations()
		#---Save the data
		pickledump(mset,'pkl.structures.'+sysname+'.'+basename[:11]+'.'+str(timeslice[0])+'-'+
			str(timeslice[1])+'-'+str(timeslice[2])+'.pkl',directory=pickles)
		if erase_when_finished:
			del mset
		print 'time = '+str(1./60*(time.time()-starttime))+' minutes.'


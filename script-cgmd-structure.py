#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = 1
framecount = None
location = 'dark'
execfile('locations.py')

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'

#---Analysis plan
analysis_descriptors = [
	(['membrane-v623-stress-test'],director_cgmd,selector_cgmd,-2),
	(['membrane-v700'],director_cgmd,selector_cgmd,-1),
	(['membrane-v614'],director_cgmd,selector_cgmd,0),
	(['membrane-v599'],director_cgmd,selector_cgmd,-1),
	(['membrane-v612'],director_cgmd,selector_cgmd,-1),
	(['membrane-v598'],director_cgmd,selector_cgmd,-1),
	(['membrane-v595'],director_cgmd,selector_cgmd,-1),
	(['membrane-v596'],director_cgmd,selector_cgmd,-1),
	(['membrane-v550'],['name PO4','name C2A'],selector_cgmd,-1),
	(['membrane-v032'],['name PO4','name C2A'],selector_cgmd,-1)]

#---Functions
#-------------------------------------------------------------------------------------------------------------

def analyze_structure(testno,traj):
	'''Compute the average structure and fluctuations of a CGMD bilayer.'''
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_cgmd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='cgmd')
	#---Average structure calculation
	mset.identify_monolayers(director,startframeno=1)
	mset.midplaner(selector,skip=skip,rounder=20.0,framecount=framecount,protein_selection='name BB')
	mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
	mset.analyze_undulations(redundant=0)
	#---Save the data
	pickledump(mset,'pkl.avgstruct.'+tests[testno]+'.'+basename+'.pkl',directory=pickles)
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors[-1:]:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,trajno) = ad
	for t in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[t]+'.'
		for traj in [trajectories[systems.index(tests[t])][trajno]]:
			#---Run the analysis function on the desired system
			mset = analyze_structure(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


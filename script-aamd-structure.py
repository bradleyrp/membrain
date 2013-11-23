#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = 1
framecount = 15
location = 'light'
execfile('locations.py')

#---Parameters
rounder = 4.0

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector_aamd_symmetric = 'name P'
selector_aamd_asymmetric = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
selector_aamd_asymmetric = '(name P and not resname CHL1)'
	
#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v509','membrane-v510','membrane-v511'],director_aamd_symmetric,selector_aamd_symmetric,-1),
	(['membrane-v510'],director_aamd_symmetric,'name P',-1),
	(['membrane-v511'],director_aamd_symmetric,'name P',-1),
	(['membrane-v530'],director_aamd_asymmetric,selector_aamd_asymmetric,slice(-1,None))]
	
#---Notes
# ensure that you can slice through all available xtcs in the analysis plan above
# set this up to run on dirac in slices

#---Functions
#-------------------------------------------------------------------------------------------------------------

def analyze_structure(testno,traj):
	'''Compute the average structure and fluctuations of an AAMD bilayer.'''
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_aamd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='cgmd')
	#---Average structure calculation
	mset.identify_monolayers(director,startframeno=1)
	mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount)
	mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
	mset.analyze_undulations(redundant=0)
	#---Save the data
	pickledump(mset,'pkl.structures.'+tests[testno]+'.'+basename+'.pkl',directory=pickles)
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors[analysis_plan]:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,trajno) = ad
	for t in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[t]+'.'
		for traj in trajectories[systems.index(tests[t])][trajno]:
			#---Run the analysis function on the desired system
			mset = analyze_structure(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


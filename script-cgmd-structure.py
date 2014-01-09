#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = None
framecount = None
location = ''
execfile('locations.py')

#---Parameters
rounder = 20.0

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[360000,460000,200]),
	(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[60000,160000,200]),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[100000,200000,200]),
	(['membrane-v550'],director_cgmd,selector_cgmd,None,slice(-1,None),[300000,400000,200])]
	
#---Functions
#-------------------------------------------------------------------------------------------------------------

def analyze_structure(testno,traj):
	'''Compute the average structure and fluctuations of a CGMD bilayer.'''
	mset = MembraneSet()
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_cgmd_surfacer
	print 'Accessing '+basename+'.'
	starttime = time.time()
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='cgmd')
	print 'Load complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'
	#---Average structure calculation
	mset.identify_monolayers(director,startframeno=0)
	if protein_select == None:
		mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,timeslice=timeslice)
	else:
		mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,protein_selection=protein_select,timeslice=timeslice)
	mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
	mset.analyze_undulations(redundant=0)
	#---Save the data
	pickledump(mset,'pkl.structures.'+tests[testno]+'.'+basename[:11]+'.'+str(timeslice[0])+'-'+
		str(timeslice[1])+'-'+str(timeslice[2])+'.pkl',directory=pickles)
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors[analysis_plan]:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,protein_select,trajno,timeslice) = ad
	for t in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[t]+'.'
		for traj in trajectories[systems.index(tests[t])][trajno]:
			#---Run the analysis function on the desired system
			mset = analyze_structure(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


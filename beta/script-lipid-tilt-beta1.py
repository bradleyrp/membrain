#!/usr/bin/python -i

from membrainrunner import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = 5
framecount = None
location = ''
execfile('locations.py')

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v599-select'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v550'],director_cgmd,selector_cgmd,None,slice(-1,None))]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def analyze_tilt(testno,traj):
	'''Compute the surface normals and lipid tilts in the top monolayer of a CGMD bilayer.'''
	mset = MembraneSet()
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_cgmd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='cgmd')
	#---Average structure calculation
	mset.identify_monolayers(director,startframeno=0)
	if protein_select == None:
		mset.tilter(selector,director,skip=skip,framecount=framecount)
	else:
		mset.tilter(selector,director,skip=skip,framecount=framecount,protein_selection=protein_select)
	#---Save the data
	pickledump(mset,'pkl.tilt.'+tests[testno]+'.'+basename+'.pkl',directory=pickles)
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors[analysis_plan]:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,protein_select,trajno) = ad
	for t in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[t]+'.'
		for traj in trajectories[systems.index(tests[t])][trajno]:
			#---Run the analysis function on the desired system
			mset = analyze_tilt(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


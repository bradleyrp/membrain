#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------


############ I htink this is deprecated. see script-dimple.py




#---Analysis parameters
location = 'light'
skip = 1
framecount = None

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-repository/'
	locations = '/home/rpb/worker/membrain/trajectory-rpb-light'	
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/store-delta/worker/worker-big/membrane-repository/'
	locations = '/store-delta/worker/membrain/trajectory-rpb-dark'
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)	

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'

#---Analysis plan
analysis_descriptors = [
	(['membrane-v700'],director_cgmd,selector_cgmd,-1)]
	
#---Analysis plan
analysis_descriptors = [
	(['membrane-v614'],director_cgmd,selector_cgmd,-1)]

#---Analysis plan
analysis_descriptors = [
	(['membrane-v612'],director_cgmd,selector_cgmd,-1)]

#---Functions
#-------------------------------------------------------------------------------------------------------------

def compute_meshes(testno,traj):
	'''Generate meshes of the top monolayer for a trajectory. Necessary for computing curvature.'''
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_cgmd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='cgmd')
	#---Average structure calculation
	mset.identify_monolayers(director,startframeno=1)
	mset.mesher(selector,skip=skip,framecount=framecount,protein_selection='name BB')
	return mset
	
def compute_curvatures(testno,traj):
	'''Compute curvatures from pre-computed meshes.'''
	#---Load the trajectory
	#---Save the data
	pickledump(mset,'pkl.mesher.'+tests[testno]+'.'+basename+'.pkl')
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,trajno) = ad
	for t in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[t]+'.'
		for traj in [trajectories[systems.index(tests[t])][trajno]]:
			#---Run the analysis function on the desired system
			mset = compute_meshes(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


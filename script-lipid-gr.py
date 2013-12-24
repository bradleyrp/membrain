#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = None
framecount = None
location = ''
execfile('locations.py')

#---Selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---Analysis plan
analysis_descriptors = [
	(['membrane-v510'],'all',director_symmetric,-1,
		[['resname DOPC and name P','name CL','DOPC P-to-CL'],
		['resname DOPC and name P','name MG','DOPC P-to-MG']]),
	(['membrane-v530'],'all',director_symmetric,-1,
		[['resname DOPS and name P','resname PI2P and name P','DOPS-PIP2']])]

#---Functions
#-------------------------------------------------------------------------------------------------------------

(tests,selector,director,trajno,pairs) = analysis_descriptors[-1]
traj = [trajectories[systems.index(tests[t])][trajno]][-1]
t = 0

#def analyze_ion_distributions(testno,traj):
if 1:
	testno = 
	'''Compute the lipid-lipid g(r).'''
	mset = MembraneSet()
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_aamd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='aamd')
	director = ['name P','name C218','name C318']
	mset.identify_monolayers(director)
	#---Batches of g(r)-voronoi calculations
	for pair in pairs:
		print pair
		mset.batch_gr_lipid_lipid([pair[0],pair[1]],framecount=framecount,skip=skip,
			label=pair[2],mode='voronoi_bin',monolayer_rep='P')
	#---Save the data
	pickledump(mset,'pkl.gr-vornoi.'+tests[testno]+'.'+basename+'.pkl')
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

'''
starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,selector,director,trajno,pairs) = ad
	for t in range(len(tests)):
		print 'Running calculation: monolayer unstructured triangulation '+tests[t]+'.'
		for traj in [trajectories[systems.index(tests[t])][trajno]]:
			#---Run the analysis function on the desired system
			mset = analyze_ion_distributions(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'
'''

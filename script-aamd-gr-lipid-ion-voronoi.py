#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
location = 'light'
skip = None
framecount = 30

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-repository'
	locations = '/home/rpb/worker/membrain/trajectory-rpb-light'	
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/'
	locations = '/store-delta/worker/membrain/trajectory-rpb-dark'
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)	

#---Selections
#---Note: this code also uses the following globals: sel_aamd_surfacer
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---Analysis plan ########### old style
analysis_descriptors = [
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-25),
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v514'],['NA'],['DOPC','DOPS','PIPU'],'all',director_symmetric,-1),
	(['membrane-v515'],['NA'],['DOPC','DOPS','PIPP'],'all',director_symmetric,-1),
	(['membrane-v510','membrane-v511'],['MG','Cal'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v530','membrane-v531','membrane-v532'],['NA','MG','Cal'],['POPC','CHL1','DOPE','DOPS','PIP2'],
		'all',director_asymmetric,-1)]
		
#---Analysis plan
analysis_descriptors = [
	(['membrane-v510'],'all',director_symmetric,-1,
		[['resname DOPC and name P','name CL','DOPC P-to-CL'],
		['resname DOPC and name P','name MG','DOPC P-to-MG']])]

#---Functions
#-------------------------------------------------------------------------------------------------------------

def analyze_ion_distributions(testno,traj):
	'''Compute the average structure and fluctuations of a CGMD bilayer.'''
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
		mset.batch_gr_lipid_ion([pair[0],pair[1]],framecount=framecount,skip=skip,
			label=pair[2],mode='voronoi_bin',monolayer_rep='P')
	#---Save the data
	pickledump(mset,'pkl.gr-vornoi.'+tests[testno]+'.'+basename+'.pkl')
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

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


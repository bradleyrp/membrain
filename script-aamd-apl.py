#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
location = 'dirac'
skip = None
framecount = 30

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-test-set'
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

#---Analysis plan
analysis_descriptors = [
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-25),
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v514'],['NA'],['DOPC','DOPS','PIPU'],'all',director_symmetric,-1),
	(['membrane-v515'],['NA'],['DOPC','DOPS','PIPP'],'all',director_symmetric,-1),
	(['membrane-v510','membrane-v511'],['MG','Cal'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v530','membrane-v531','membrane-v532'],['NA','MG','Cal'],['POPC','CHL1','DOPE','DOPS','PIP2'],
		'all',director_asymmetric,-1)]

#---Functions
#-------------------------------------------------------------------------------------------------------------

def analyze_area_per_lipid(testno,traj):
	'''For a particular trajectory, calculate the Voronoi cells and areas-per-lipid.'''
	mset = MembraneSet()
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	print 'Accessing '+basename+'.'
	mset.load_trajectory(
		sel_aamd_surfacer,
		(basedir+'/'+gro,basedir+'/'+traj),
		resolution='aamd')
	mset.identify_monolayers(director,startframeno=3)
	mset.identify_residues(residues)
	#---Area-per-lipid calculations
	mset.triangulator(selector,framecount=framecount,skip=skip,tesstype='voronoi')
	#---Save the data
	pickledump(mset,'pkl.apl-cells.'+tests[testno]+'.'+basename+'.pkl')
	#---Clear the class object
	return mset	

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,ionnames,residues,selector,director,trajno) = ad
	for t in range(len(tests)):
		print 'Running calculation: monolayer unstructured triangulation '+tests[t]+'.'
		for traj in [trajectories[systems.index(tests[t])][trajno]]:
			#---Run the analysis function on the desired system
			mset = analyze_area_per_lipid(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


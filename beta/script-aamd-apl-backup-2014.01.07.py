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
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-25),
	(['membrane-v509'],['NA'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v514'],['NA'],['DOPC','DOPS','PIPU'],'all',director_symmetric,-1),
	(['membrane-v515'],['NA'],['DOPC','DOPS','PIPP'],'all',director_symmetric,-1),
	(['membrane-v510','membrane-v511'],['MG','Cal'],['DOPC','DOPS','PI2P'],'all',director_symmetric,-1),
	(['membrane-v530','membrane-v531','membrane-v532'],['NA','MG','Cal'],['POPC','CHL1','DOPE','DOPS','PIP2'],
		'all',director_asymmetric,-1),
	(['membrane-v534'],['Cal'],['POPC','CHL1','DOPE','DOPS','P35P'],
		'all',director_asymmetric,slice(-1,None)),
	(['membrane-v533'],['Mg'],['POPC','CHL1','DOPE','DOPS','P35P'],
		'all',director_asymmetric,slice(-1,None))]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors[analysis_plan]:
	(tests,ionnames,residues,selector,director,trajno) = ad
	for testno in range(len(tests)):
		print 'Running calculation: monolayer unstructured triangulation '+tests[testno]+'.'
		for traj in trajectories[systems.index(tests[testno])][trajno]:
			#---Area-per-lipid analysis routing (note that this is a membrain.py built-in.
			mset = MembraneSet()
			#---Load the trajectory
			gro = structures[systems.index(tests[testno])]
			basename = traj.split('/')[-1][:-4]
			print 'Accessing '+basename+'.'
			mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director,startframeno=0)
			mset.identify_residues(residues)
			#---Area-per-lipid calculations
			mset.triangulator(selector,framecount=framecount,skip=skip,tesstype='voronoi')
			#---Save the data
			pickledump(mset,'pkl.apl-cells.'+tests[testno][9:14]+'.'+basename+'.pkl',directory=pickles)
			#---Clear the class object
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


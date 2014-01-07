#!/usr/bin/python

#---Settings
#-------------------------------------------------------------------------------------------------------------

from membrainrunner import *

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
analysis_plan = slice(-2,None)
analysis_descriptors = [
	(['membrane-v533'],'all',director_asymmetric,slice(-2,-1),
		['POPC','CHL1','DOPE','DOPS','P35P'],
		[['resname P35P and name P','name CL','P35P P-to-CL'],
		['resname P35P and name P','name MG','P35P P-to-MG']],'P35P-to-ions'),
	(['membrane-v534'],'all',director_asymmetric,slice(-2,-1),
		['POPC','CHL1','DOPE','DOPS','P35P'],
		[['resname P35P and name P','name CL','P35P P-to-CL'],
		['resname P35P and name P','name Cal','P35P P-to-CA']],'P35P-to-ions'),
	(['membrane-v531'],'all',director_asymmetric,slice(-2,-1),
		['POPC','CHL1','DOPE','DOPS','PI2P'],
		[['resname PI2P and name P','name CL','PI2P P-to-CL'],
		['resname PI2P and name P','name MG','PI2P P-to-MG']],'PI2P-to-ions'),
	(['membrane-v532'],'all',director_asymmetric,slice(-2,-1),
		['POPC','CHL1','DOPE','DOPS','PI2P'],
		[['resname PI2P and name P','name CL','PI2P P-to-CL'],
		['resname PI2P and name P','name Cal','PI2P P-to-CA']],'PI2P-to-ions')]

#---Methods
do_compute_lipid_ion_distribution_cells = 1

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Compute lipid-ion distributions via Vornoi cells method
if do_compute_lipid_ion_distribution_cells:
	starttime = time.time()
	#---loop over analysis descriptors
	for ad in analysis_descriptors[analysis_plan]:
		(tests,selector,director,trajno,residues,pairs,extraname) = ad
		#---loop over tests within the descriptor
		for testno in range(len(tests)):
			#---loop over specified trajectories
			for traj in trajectories[systems.index(tests[testno])][trajno]:
				mset = MembraneSet()
				gro = structures[systems.index(tests[testno])]
				basename = traj.split('/')[-1][:-4]
				print 'Accessing '+basename+'.'
				mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
				mset.identify_monolayers(director,startframeno=0)
				mset.identify_residues(residues)
				#---frame selection header
				end = None
				start = None
				if framecount == None:
					if end == None: end = mset.nframes
					if start == None: start = 0
					if skip == None: skip = 1
				else:
					start = 0
					end = mset.nframes
					skip = int(float(mset.nframes)/framecount)
					skip = 1 if skip < 1 else skip
				for pair in pairs:
					print pair
					mset.batch_gr_lipid_ion([pair[0],pair[1]],framecount=framecount,skip=skip,start=start,
						end=end,label=pair[2],mode='voronoi_bin',monolayer_rep='P')
				#---Save the data
				pickledump(mset,'pkl.gr-vornoi.'+tests[testno][9:14]+'.'+basename+'.'+
					extraname+'.pkl',directory=pickles)


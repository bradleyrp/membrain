#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')

#---selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector_aamd_symmetric = 'name P'
selector_aamd_asymmetric = '(name P and not resname CHL1)'
selector_aamd_asymmetric = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
residues_aamd_symmetric = ['DOPC','DOPS','PI2P']
residues_aamd_asymmetric = ['DOPC','DOPS','DOPE','POPC','PI2P', 'CHL1', 'P35P']
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(None,None)
analysis_descriptors = {
	'v533-40000-54000-100':
		{'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-kraken-md.part0015.40000-54000-100.atomP.xtc',
		'whichframes':None},
	'v534-40000-60000-100':
                {'sysname':'membrane-v534',
                'sysname_lookup':'membrane-v534-atomP',
                'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
                'residues':'infer',
                'trajsel':'s4-sim-kraken-md.part0013.40000-60000-100.atomP.xtc',
                'whichframes':None},
	'v530-30000-100000-100':
		{'sysname':'membrane-v530',
                'sysname_lookup':'membrane-v530-atomP',
                'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
                'residues':'infer',
                'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.atomP.xtc',
                'whichframes':None},
	'v531-20000-62000-100':
		{'sysname':'membrane-v531',
                'sysname_lookup':'membrane-v531-atomP',
                'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
                'residues':'infer',
                'trajsel':'s4-sim-trestles-md.part0007.20000-62000-100.atomP.xtc',
                'whichframes':None},
	'v532-20000-58000-100':
		{'sysname':'membrane-v532',
                'sysname_lookup':'membrane-v532-atomP',
                'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
                'residues':'infer',
                'trajsel':'s4-sim-trestles-md.part0007.20000-58000-100.atomP.xtc',
                'whichframes':None},
	}
analysis_names = [
	'v530-30000-100000-100',
	'v531-20000-62000-100',
	'v532-20000-58000-100'
	][:]

#---MAIN
#-------------------------------------------------------------------------------------------------------------
#---loop over analysis questions
for aname in analysis_names:
	#---details
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
	#---loop over trajectory files
	for traj in trajfile:
		mset = MembraneSet()
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		mset.identify_monolayers(director)
		if residues == 'infer':
			residues = list(set(mset.universe.residues.resnames()))
		mset.identify_residues(residues)
		result_data = MembraneData('cells')
		if whichframes == None: whichframes = slice(None,None)
		#---calcalate areas per lipid
		for frameno in range(mset.nframes)[whichframes]:
			print frameno
			#---select the frame
			mset.gotoframe(frameno)
			if type(selector) == str:
				selector = [selector for i in range(len(mset.resnames))]
			#---collect the lipid positions
			top = []
			for t in range(len(mset.resnames)):
				top.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
					axis=0) for i in mset.monolayer_by_resid[0][t]])
			bot = []
			for t in range(len(mset.resnames)):
				bot.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
					axis=0) for i in mset.monolayer_by_resid[1][t]])
			#---calculate areas for each monolayer
			results = []
			for points in [top,bot]:
				points_pbc = mset.wrappbc(points,vecs=mset.vec(frameno),mode='grow')
				vmap = scipy.spatial.Voronoi(points_pbc[:,0:2])
				areas = []
				for p in range(len(points)):
					vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
					pairs = zip(vertices, vertices[1:]+vertices[0:1])
					areas.append(abs(sum(x1*y2-y1*x2 for (x1, y1), (x2, y2) in pairs)/2))
				results.append([points,vmap,areas])
			result_data.add(results,[frameno])
			#---add details to the store
			for i in analysis_descriptors[aname]:
				result_data.addnote([i,(analysis_descriptors[aname])[i]])
			result_data.addnote(['selector',selector])
			result_data.addnote(['director',director])
		mset.store.append(result_data)
		#---Save the data
		pickledump(mset,'pkl.cells.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)

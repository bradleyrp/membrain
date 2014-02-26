#!/usr/bin/python -i

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
analysis_plan = slice(None,None)
analysis_descriptors = {
	'v530-5000-5500-2':
		{'sysname':'membrane-v530-apl',
		'sysname_lookup':'membrane-v530',
		'director':director_aamd_symmetric,'selector':selector_aamd_symmetric,'protein_select':None,
		'trajsel':'u4-reboot-sim-compbio-md.part0002.5000-5500-2.xtc'},
	'v530-20000-20500-2':
		{'sysname':'membrane-v530-apl',
		'sysname_lookup':'membrane-v530',
		'director':director_aamd_symmetric,'selector':selector_aamd_symmetric,'protein_select':None,
		'trajsel':'u4-reboot-sim-compbio-md.part0005.20000-20500-2.xtc'},
	'v530-80000-80500-2':
		{'sysname':'membrane-v530-apl',
		'sysname_lookup':'membrane-v530',
		'director':director_aamd_symmetric,'selector':selector_aamd_symmetric,'protein_select':None,
		'trajsel':'u4-reboot-sim-compbio-md.part0005.80000-80500-2.xtc'}
		
		}
analysis_names = ['v530-5000-5500-2', 'v530-20000-20500-2', 'v530-80000-80500-2']


#---Methods
do_compute_area_cells = 1

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Compute lipid areas via Vornoi cells method
if do_compute_area_cells:
	starttime = time.time()
	#---loop over analysis descriptors
	for ad in analysis_descriptors[analysis_plan]:
		(tests,ionnames,residues,selector,director,trajno) = ad
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
				result_data = MembraneData('cells')
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
				#---calcalate areas per lipid
				for frameno in range(start,end,skip):
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
				mset.store.append(result_data)
				#---Save the data
				pickledump(mset,'pkl.apl-cells.'+tests[testno][9:14]+'.'+basename+'.pkl',directory=pickles)


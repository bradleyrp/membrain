#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

if 'batch_override' not in globals():

	execfile('header-cgmd.py')
	analysis_names = [
		'v620-50000-80000-50',
		][-1:]
	routine = [
		'calc_apl',
		'plot_apl',
		'mov',
		][:1]
	flaglist = []
	
#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'calc_apl' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			#---note: removed a check to see if the pkl exists from this space
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			result_data_points = MembraneData('cells')
			result_data_vmap = MembraneData('cells')
			result_data_areas = MembraneData('cells')
			if whichframes == None: whichframes = slice(None,None)
			#---calcalate areas per lipid
			for frameno in range(mset.nframes)[whichframes]:
				status('status: storing points for fr = ',i=frameno,looplen=mset.nframes)
				#---select the frame
				mset.gotoframe(frameno)
				if type(selector) == str:
					selector = [selector for i in range(len(mset.resnames))]
				#---collect the lipid positions
				top = []
				for t in range(len(mset.resnames)):
					top.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid_abs[0][t]])
				bot = []
				for t in range(len(mset.resnames)):
					bot.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid[1][t]])
				#---calculate areas for each monolayer
				results_points = []
				results_vmap = []
				results_areas = []
				for points in [top,bot]:
					points_pbc = mset.wrappbc(points,vecs=mset.vec(frameno),mode='grow')
					vmap = scipy.spatial.Voronoi(points_pbc[:,0:2])
					areas = []
					for p in range(len(points)):
						vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
						pairs = zip(vertices, vertices[1:]+vertices[0:1])
						areas.append(abs(sum(x1*y2-y1*x2 for (x1, y1), (x2, y2) in pairs)/2))
					results_points.append([points,[],[]])
					results_vmap.append([[],vmap,[]])
					results_areas.append([[],[],areas])
				result_data_points.add(results_points,[frameno])
				result_data_vmap.add(results_vmap,[frameno])
				result_data_areas.add(results_areas,[frameno])
			#---add details to the store
			for result_data in [result_data_points,result_data_vmap,result_data_areas]:
				for i in analysis_descriptors[aname]:
					result_data.addnote([i,(analysis_descriptors[aname])[i]])
				result_data.addnote(['selector',selector])
				result_data.addnote(['director',director])
			results_table = [
				['points',result_data_points],
				['vmap',result_data_vmap],
				['areas',result_data_areas],
				]
			for thisres in results_table:
				mset.store = [thisres[1]]
				#---Save the data
				pickledump(mset,'pkl.cells.'+thisres[0]+'-'+'-'.join(flaglist)+('-' if len(flaglist) > 0 else '')+specname_pickle(sysname,traj)+'.pkl',directory=pickles)
			checktime()

if 'mov' in routine:

	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in range(9)]
	clrset = [clrs[i] for i in [8,7,1,0]]

	#---relative resids
	resids = [[mset.resids_reracker.index(i) for i in j] for j in mset.resids]
	#---color codes
	colorcodes = [[[i for i in range(4) if j in resids[i]][0] 
		for j in mono] for mono in mset.monolayer_residues]

	#---populate the data structure for the movie maker
	dat_vor = [[0. for cols in range(1)] for rows in range(2)]
	#---rows are monolayers, columns are systems
	dat_vor[0][0] = list(dat.get(['monolayer',0,'type','voronoi']))
	dat_vor[1][0] = list(dat.get(['monolayer',1,'type','voronoi']))
	rand_color_list = [np.random.rand(3,1) for i in range(2*len(dat_vor[0][0][0].points))]
	#---celltest was updated slightly, see the plot_snapshot section
	plotmov(dat_vor,'celltest',plotfunc='cellplot')
	

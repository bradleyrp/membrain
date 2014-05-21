#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')

execfile('header-aamd.py')

#---settings
compare_phosphate_position = False

#---standard selection
analysis_names = [
	'v530-40000-90000-50',
	'v531-40000-90000-50',
	'v532-40000-90000-50',
	'v509-40000-90000-50',
	'v510-40000-90000-50',
	'v511-40000-90000-50',
	'v533-40000-90000-50',
	'v534-40000-90000-50',
	][:]
	
#---alternate tests
if compare_phosphate_position:
	analysis_names = [
		'v531-40000-90000-50',
		'v532-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		]

routine = [
	'calc_apl',
	'plot_apl',
	'calc_span',
	'plot_span',
	][-2:-1]
	
#---selector overrides
selector_override = [
	'no_chl',
	None,
	][0]
	
#---which types of plots to generate
resname_group = [
	'all',
	'phosphate_position',
	][0]
	
if compare_phosphate_position: resname_group = 'phosphate_position'

#---custom settings
flaglist = []
if selector_override == 'no_chl':
	for aname in analysis_names:
		(analysis_descriptors[aname])['residues'] = residues_aamd_asymmetric_no_chl	
		(analysis_descriptors[aname])['director'] = director_aamd_asymmetric_no_chl
		(analysis_descriptors[aname])['selector'] = selector_aamd_asymmetric_no_chl
	flaglist.append('no_chl')

#---plot settings
nbins = 40
showplots = False
hist_area_max = 140

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
				status('status: storing points for fr = '+str(frameno+1)+'/'+str(mset.nframes))
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
				pickledump(mset,'pkl.cells.'+thisres[0]+'-'+'-'.join(flaglist)+\
					('-' if len(flaglist) > 0 else '')+\
					specname_pickle(sysname,traj)+'.pkl',directory=pickles)
			checktime()

if 'calc_span' in routine:
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
			#---empties
			head_area = []
			head_angle = []
			data_max = []
			mdat = []
			result_data = []
			result_data = MembraneData('spanangle2')
			residues = mset.universe.selectAtoms('resname '+ptdins_resname).resids()
			whichframes = range(len(mset.universe.trajectory))
			for fr in whichframes:
				mset.gotoframe(fr)		
				status('status: frame = '+str(fr+1)+'/'+str(len(whichframes)))
				results_by_frame = []
				for res in residues:
					coords = mset.universe.selectAtoms(\
						'resid '+str(res)+' and resname '+ptdins_resname+' and '+headspan)
					head_size = (max(scipy.spatial.distance.pdist(coords.coordinates())))
					head_area.append(head_size*head_size)
					coords = mset.universe.selectAtoms(\
						'resid '+str(res)+' and '+headangle).coordinates()
					angle = (arccos(np.dot(coords[0]-coords[1],coords[2]-coords[1])/\
						np.linalg.norm(coords[0]-coords[1])/np.linalg.norm(coords[2]-coords[1])))
					head_angle.append(angle*(180./pi))
					results_by_frame.append([
						head_size*head_size,
						angle*(180./pi),
						mset.universe.trajectory[fr].time,
						str(res)
						])
				result_data.add(results_by_frame,[fr])
			#---descriptors in the pkl
			for i in analysis_descriptors[aname]:
				result_data.addnote([i,(analysis_descriptors[aname])[i]])
			#---save	
			mset.store.append(result_data)
			pickle.dump(mset,open(pickles+'pkl.headspan-headangle2.'+\
				specname_pickle(sysname,traj)+'.pkl','w'))
			checktime()

#---LOADS
#-------------------------------------------------------------------------------------------------------------

if 'calc_apl' not in routine and 'plot_apl' in routine and 'msets' not in globals():
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.cells.'+\
			'areas-'+'-'.join(flaglist)+('-' if len(flaglist) > 0 else '')+\
			specname_pickle(sysname,trajfile[0])+'.pkl')
		msets.append(mset)
		
if 'calc_span' not in routine and 'plot_span' in routine and 'msets' not in globals():
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.headspan-headangle2.'+aname+'.pkl')
		msets.append(mset)

#---PLOT
#-------------------------------------------------------------------------------------------------------------

if 'plot_apl' in routine:

	#---two-dimensional list of which residues to compare across simulations
	if resname_group == 'all':
		plot_resnames = list(set([j for k in [i.resnames for i in msets] for j in k]))
	elif resname_group == 'phosphate_position':
		plot_resnames = ['ptdins']

	#---plot on a single panel
	largest_area = max([
		max(flatten(
		msets[analysis_names.index(aname)].getdata('cells').get(['type','areas'])
		)) for aname in analysis_names])
	#---note: largest areas was giving one that was 1400 A2, clearly an error which needs checked
	largest_area = hist_area_max

	#---loop over residue comparisons
	for resname in plot_resnames:
		lims = (0,largest_area)
		fig = plt.figure()
		ax = plt.subplot(111)
		peakval = 0
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			mset = msets[analysis_names.index(aname)]
			if plot_resnames == ['ptdins']: resname = ptdins_resname
			if resname in mset.resnames:
				data = []
				for mono in range(2):
					inds = array([mset.monolayer_residues[mono].index(i) 
						for i in mset.monolayer_by_resid[mono][mset.resnames.index(resname)]])
					if len(inds) > 0:
						data.extend(\
							[array(i)[inds] for i in \
							mset.getdata('cells').get(['monolayer',mono,'type','areas'])]
							)
				hist,binedges = histogram(data,bins=nbins,normed=True,range=lims)
				mids = (binedges[1:]+binedges[:-1])/2
				peakval = max(hist) if max(hist) > peakval else peakval
				if resname == ptdins_resname: label = ptdins_label
				else: label = resname
				if resname_group == 'phosphate_position':
					color = color_dictionary_aamd(ionname=ion_name,lipid_resname=resname,
						comparison='ions_phospate_position')
				else:
					color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
				#---plot commands
				ax.plot(mids,hist,'-',lw=2,
					color=color)
				fill_between(mids,[0 for i in range(len(hist))],hist,
					facecolor=color,
					alpha=0.35,
					label=label+', '+ion_label,
					ax=ax)
		#---plot settings
		ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
		ax.grid(True)
		ax.legend()
		ax.set_ylim((0,peakval*1.2))
		plt.setp(ax.get_xticklabels(),fontsize=fsaxticks)
		ax.set_ylabel('relative frequency',fontsize=fsaxlabel)
		ax.set_xlabel(r'area per molecule $(\mathrm{\AA^2})$',fontsize=fsaxlabel)
		#---save
		plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-apl-'+resname+'-'+\
			('no_chl-' if selector_override == 'no_chl' else '')+\
			'-'.join(analysis_names)+'.png',\
			dpi=300,bbox_inches='tight')
		if showplots: plt.show()
		plt.close(fig)

if 'plot_span' in routine:

	largest_area = hist_area_max
	lims = (0,largest_area)
	fig = plt.figure()
	ax = plt.subplot(111)
	peakval = 0
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		mset = msets[analysis_names.index(aname)]
		tmp = mset.getdata('spanangle2')
		tmp2 = tmp.get(['type','headspan'])
		#---note due to a slight issue with the membraindata get function, this returns array of strings
		data = [float(i) for j in tmp2 for i in j]
		hist,binedges = histogram(data,bins=nbins,normed=True,range=lims)
		mids = (binedges[1:]+binedges[:-1])/2
		peakval = max(hist) if max(hist) > peakval else peakval
		if resname_group == 'phosphate_position':
			color = color_dictionary_aamd(ionname=ion_name,lipid_resname=ptdins_resname,
				comparison='ions_phospate_position')
		else:
			color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
		#---plot commands
		ax.plot(mids,hist,'-',lw=2,
			color=color)
		fill_between(mids,[0 for i in range(len(hist))],hist,
			facecolor=color,
			alpha=0.35,
			label=ptdins_label+', '+ion_label,
			ax=ax)

	#---plot settings
	ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
	ax.grid(True)
	ax.legend(ncol=2)
	ax.set_ylim((0,peakval*1.3))
	plt.setp(ax.get_xticklabels(),fontsize=fsaxticks)
	ax.set_ylabel('relative frequency',fontsize=fsaxlabel)
	ax.set_xlabel(r'area per molecule $(\mathrm{\AA^2})$',fontsize=fsaxlabel)
	#---save
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-span-'+\
		'-'.join(analysis_names)+'.png',\
		dpi=300,bbox_inches='tight')
	if showplots: plt.show()
	plt.close(fig)
	
#---DEVELOPMENT
#-------------------------------------------------------------------------------------------------------------

def cellplot(dat,fig,altdat=None,vmin=None,vmax=None,fr=None,panels=1):
	rows,cols = shape(dat)
	gs = gridspec.GridSpec(rows,cols,wspace=0.0,hspace=0.0)
	#---rows are monolayers, columns are systems
	for row in range(rows):
		for col in range(cols):
			vor = dat[row][col]
			ax = fig.add_subplot(gs[row,col],aspect='equal')
			regionlist = [vor.regions[i] for i in vor.point_region]
			for r in range(len(regionlist)):
				region = regionlist[r]
				if not -1 in region:
					polygon = [vor.vertices[i] for i in region]
					if 0: axes.fill(*zip(*polygon),alpha=0.5)
					p = mpl.patches.Polygon(polygon,alpha=0.65,
						facecolor=('w' if r >= len(colorcodes[row]) else clrset[colorcodes[row][r]]),
						lw=0.5,edgecolor='k')
					ax.add_patch(p)
			for p in range(len(vor.points)):
				pt = vor.points[p]
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1.,facecolor=rand_color_list[p],
					edgecolor='k',lw=0.3))
			ax.set_xlim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][0])
			ax.set_ylim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][1])
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xticks([])
			ax.set_yticks([])
	return [gs,ax]

#---code for making grid movie, under construction
if 'plot_mov' in routine:

	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in range(9)]
	clrset = [clrs[i] for i in [8,7,1,0]]

	def cellplot(dat,fig,altdat=None,vmin=None,vmax=None,fr=None,panels=1):
		rows,cols = shape(dat)
		gs = gridspec.GridSpec(rows,cols,wspace=0.0,hspace=0.0)
		#---rows are monolayers, columns are systems
		for row in range(rows):
			for col in range(cols):
				vor = dat[row][col]
				ax = fig.add_subplot(gs[row,col],aspect='equal')
				regionlist = [vor.regions[i] for i in vor.point_region]
				for r in range(len(regionlist)):
					region = regionlist[r]
					if not -1 in region:
						polygon = [vor.vertices[i] for i in region]
						if 0: axes.fill(*zip(*polygon),alpha=0.5)
						p = mpl.patches.Polygon(polygon,alpha=0.65,
							facecolor=('w' if r >= len(colorcodes[row]) else clrset[colorcodes[row][r]]),
							lw=0.5,edgecolor='k')
						ax.add_patch(p)
				for p in range(len(vor.points)):
					pt = vor.points[p]
					ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1.,facecolor=rand_color_list[p],
						edgecolor='k',lw=0.3))
				ax.set_xlim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][0])
				ax.set_ylim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][1])
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				ax.set_xticks([])
				ax.set_yticks([])
		return [gs,ax]

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
	#plotmov(dat_vor,'celltest',plotfunc='cellplot')




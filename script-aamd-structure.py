#!/usr/bin/python

logfile,interact,debugmode = [None,False,None]
from membrainrunner import *
execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---parameters
rounder = 4

#---load the standard header defintions
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
	][:3]
	
#---alternate tests
if compare_phosphate_position:
	analysis_names = [
		'v531-40000-90000-50',
		'v532-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		]
	
#---routine
routine = [
	'calc',
	'load',
	'plot',
	][2:]
	
#---settings
showplots = False

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'calc' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			isfile = unpickle(pickles+'pkl.structures.'+'space'+str(rounder)+'A.'+\
				specname_pickle(sysname,traj)+'.pkl')
			if isfile == True:
				status('status: pkl already exists ')
			else:
				mset = MembraneSet()
				#---load the trajectory
				status('status: accessing '+specname_pickle(sysname,traj))
				mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
				checktime()
				#---average structure calculation
				mset.identify_monolayers(director)
				#---infer the timeslice from the XTC name if not specified
				if 'timeslice' in analysis_descriptors[aname].keys():
					status('warning: requested timeslice '+str(timeslice))
				else:
					if len(trajsel.split('.')[-2].split('-')) == 1:
						tslicepos = -3
						subset_name = trajsel.split('.')[-2]
					else: tslicepos = -2
					timeslice = [int(i) for i in trajsel.split('.')[tslicepos].split('-')]
				if protein_select == None:
					mset.midplaner(selector,skip=skip,rounder=float(rounder),
						framecount=framecount,timeslice=timeslice,
						thick=True)
				else:
					mset.midplaner(selector,skip=skip,rounder=float(rounder),framecount=framecount,
						protein_selection=protein_select,timeslice=timeslice)
				mset.calculate_undulations()
				#---save the data
				pickledump(mset,'pkl.structures.'+'space'+str(rounder)+'A.'+\
					specname_pickle(sysname,traj)+'.pkl',directory=pickles)
				if erase_when_finished:
					del mset
				checktime()

#---LOAD
#-------------------------------------------------------------------------------------------------------------

if 'calc' not in routine and 'msets' not in globals():
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.structures.'+'space'+str(rounder)+'A.'+\
			specname_pickle(sysname,trajfile[0])+'.pkl')
		msets.append(mset)

#---PLOT
#-------------------------------------------------------------------------------------------------------------

if 'plot' in routine:
	fig = plt.figure()
	gs = gridspec.GridSpec(1,1)
	ax = plt.subplot(gs[0])
	for m in range(len(msets)):
		aname = analysis_names[m]
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		mset = msets[m]
		mset.calculate_undulations(peri=True)
		resname = ptdins_resname
		color = color_dictionary_aamd(ionname=ion_name,lipid_resname=resname,
			comparison='ions')
		if compare_phosphate_position:
			color = color_dictionary_aamd(ionname=ion_name,lipid_resname=resname,
				comparison='ions_phospate_position')
		plotter_undulate(mset,qmagfilter=[10**-10,1],ax=ax,inset2d=False,
			colorspec=color,showkappa=False,
			label=ptdins_label+', '+ion_label)
		if 0: plotter_undulate(mset,qmagfilter=[10**-10,1],ax=ax,inset2d=False,
			colorspec=color,showkappa=False,
			label=ptdins_label+', '+ion_label,
			peri=True)
	ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
	ax.legend(loc='upper right')
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-bending-'+\
		'-'.join(analysis_names)+'.png',\
		dpi=300,bbox_inches='tight')
	if showplots: plt.show()
	plt.close(fig)

if 0:
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
	colorcodes = [[[i for i in range(4) if j in resids[i]][0] for j in mono] for mono in mset.monolayer_residues]

	#---populate the data structure for the movie maker
	dat_vor = [[0. for cols in range(1)] for rows in range(2)]
	#---rows are monolayers, columns are systems
	dat_vor[0][0] = list(dat.get(['monolayer',0,'type','voronoi']))
	dat_vor[1][0] = list(dat.get(['monolayer',1,'type','voronoi']))
	rand_color_list = [np.random.rand(3,1) for i in range(2*len(dat_vor[0][0][0].points))]

	#plotmov(dat_vor,'celltest',plotfunc='cellplot')




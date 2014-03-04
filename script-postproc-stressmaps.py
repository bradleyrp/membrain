#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

import scipy.interpolate
import scipy.integrate
from scipy import ndimage

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Key parameters
'''
Notes:

previously the parameters were set as follows
	framewise_test = [4,32,64,1]
	test = framewise_test
	span,nnum,numgridpts,distance_factor = test
	distance_factor was used in my old pseudo-RBF function
	nnum was the number of nearest neighbors that contributed to the sum smoothed out in the pseudo-RBF
	numgridpts is the 1D length of grid points sampled from the pseudo RBF
	
other than the span, which set the number of voxels integrated on both sides of the bilayer 
and the choice of how to define the midplane
the other parameters are all smoothing parameters
hence the way to calculate the c0 maps is to do the integration, dump the then-much-smaller data
and smooth it seperately
'''

#---analysis plan
analysis_descriptors = {
	'v614-120000-220000-200': 
		{'sysname':'membrane-v614',
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'nprots':4,
		'datdir3dpp':
			'/home/rpb/compbio/membrane-v614-enthx4-12800/a8-stress-s9-120000-220000/results',
		'framewise_part':4,	
		'voxelsize':1.0,
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$'},
	'v701-60000-160000-200':
		{'sysname':'membrane-v701',
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc',
		'nprots':2,
		'datdir3dpp':
			'/home/rpb/compbio-alt/membrane-v701-exo70-anti-dilute/a1-stress-1.0-framewise/results',
		'framewise_part':3,
		'voxelsize':1.0,
		'label':r'$\textbf{{EXO70}\ensuremath{\times}2{\small (anti)}}$',
		'structure_pkl':'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
		},
	'v550-300000-400000-200': 
		{'sysname':'membrane-v550',
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'nprots':0,
		'datdir3dpp':
			'/home/rpb/compbio-alt/membrane-v550/'+\
			'a1-stress-1.0-framewise-md.part0006.300000-400000-200/results',
		'framewise_part':6,
		'voxelsize':1.0,
		'label':r'$control$'}}
analysis_names = ['v614-120000-220000-200','v550-300000-400000-200','v701-60000-160000-200'][:]
routine = ['calc_c0maps','plot','video'][1:2]
span_sweep = [1,2,3,4,5,6]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def stressmap_panel_plot(dat,fig,fr=None,cmap=None,vmax=None,vmin=None,altdat=None):
	'''Function which plots a stressmap with proteins.'''
	#---settings
	panels = len(dat)
	if cmap == None: cmap = mpl.cm.RdBu_r
	if fr == None: fr = 0
	#---axes
	gs = gridspec.GridSpec(1,panels)
	for p in range(panels):
		ax = fig.add_subplot(gs[p])
		im = ax.imshow((dat[p]).T,interpolation='nearest',origin='lower',
			cmap=cmap,vmax=vmax,vmin=vmin)
		if altdat != None and altdat[p].protein != []:
			plothull(ax,msets[p].protein[fr],griddims=shape(md_maps[0][0].data)[1:],
				vecs=mean(msets[p].vecs,axis=0),subdivide=nnprots[p],alpha=0.35,c='k')
		ax.set_title(labels[p],fontsize=fsaxtitle)
		ax.set_xlabel(r'$x\:(\mathrm{nm})$')
		ax.set_ylabel(r'$y\:(\mathrm{nm})$')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
	axins.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',
		fontsize=fsaxlabel,rotation=270)

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'calc_c0maps' in routine or ('md_maps' not in globals() and 'calc_c0maps' not in routine):
	md_maps = []
	msets = []
	labels = []
	nnprots = []
	for aname in analysis_names:
		#---get names and load the structure
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		if 'structure_pkl' not in globals() or structure_pkl == None:
			msetfile = 'pkl.structures.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		else: msetfile = structure_pkl
		mset = unpickle(pickles+msetfile)
		msets.append(mset)
		picklename = 'pkl.stressmaps.'+specname_pickle(sysname,trajfile[0])+'pkl'
		md_map = unpickle(pickles+picklename)
		labels.append(label)
		nnprots.append(nprots)
		if md_maps == None:
			print 'status: no stressmaps available so will try to compute them now'
			result_data_spans = [MembraneData('collect_c0maps') for i in range(len(span_sweep))]
			for frame in range(len(mset.surf)):
				print 'status: computing curvature maps frame: '+str(frame)
				#---load raw stress tensor data
				file3dpp = pickles+'/'+datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+\
					str('%04d'%frame)+'.lp.dat3d'
				file3dpp = datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+\
					'.fr'+str('%04d'%frame)+'.lp.dat3d'
				dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
				griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
				vecs = mset.vecs[frame]
				blocksize = (vecs/griddims)
				#---allocate space for the height-wise voxel tensions
				xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
					griddims[1]+2)])
				#---reference midplane surface
				ref_surf = mset.unzipgrid(mset.surf[frame],vecs=vecs)
				#---interpolate the reference surface to match the spacing of the stress tensors
				interp = scipy.interpolate.LinearNDInterpolator(ref_surf[:,0:2],ref_surf[:,2],fill_value=0.0)
				ref_surf_interp = array([[round((interp(i,j)+mset.surf_position[frame])/blocksize[2]) 
					for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
				#---loop over limits of integration
				for span in span_sweep:
					rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
					for pt in dat3dpp:
						if ((abs(pt[2]-ref_surf_interp[pt[0]][pt[1]]) <= span) and 
							(abs(pt[2]-ref_surf_interp[pt[0]][pt[1]]) <= span)):
							rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
								(pt[2]-ref_surf_interp[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
					#---integrate the collected stresses to determine the spontaneous curvature
					results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)]
						for i in range(griddims[0]+1)])
					#---set aside the data
					result_data_spans[span_sweep.index(span)].add(results,[frame])
			#---reorganize the data
			for r in range(len(result_data_spans)):
				result_data_spans[r].addnote(['c0_span',span_sweep[r]])
				for key in analysis_descriptors[aname]: 
					result_data_spans[r].addnote([key,analysis_descriptors[aname][key]])
				mset.store.append(result_data_spans[r])
			#---write the store, since saving mset is very costly for some reason
			#---Nb saving mset was 100MB with only 10 frames even after deleting the redundant data
			pickledump(mset.store,picklename,directory=pickles)
			md_maps = mset.store
		md_maps.append(md_map)
		

#---a single plot
if 'plot' in routine:
	#---average C0 map with PBC Gaussian blur
	if 1:
		m,n = shape(md_maps[0][0].data)[1:]
		panelplots = [scipy.ndimage.filters.gaussian_filter(numpy.tile(mean(md_maps[i][2].data,axis=0),
			(3,3)),2)[m:2*m,n:2*n] for i in range(len(md_maps))]
		vmin = array(panelplots).min()
		vmax = array(panelplots).max()
		extrem = max(abs(vmax),abs(vmin))
		vmax = extrem
		vmin = -1*extrem
		fig = plt.figure()
		stressmap_panel_plot(panelplots,fig,vmax=vmax,vmin=vmin,altdat=msets)
		plt.show()

#---make cool videos
if 'video' in routine:
	vidpanels = []
	for p in range(len(md_maps)):
		vidpanels.append(array([scipy.ndimage.filters.gaussian_filter(mean(md_maps[p][2].data[i:i+20],axis=0),2) 
			for i in range(500-20)]))
	plotmov(vidpanels,'stressmap-v614.v550.v701',altdat=msets,panels=len(md_maps),
		plotfunc='stressmap_panel_plot',whitezero=True)
		

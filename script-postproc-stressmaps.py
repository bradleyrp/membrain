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
	'v614-40000-140000-200': 
		{'sysname':'membrane-v614',
		'trajsel':'s6-sim-lonestar/md.part0002.40000-140000-200.xtc',
		'nprots':4,
		'datdir3dpp':
			'/home/rpb/compbio/membrane-v614-enthx4-12800/a9-stress-s6-40000-140000-200/results',
		'framewise_part':2,	
		'voxelsize':1.0,
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4\,v1}$'},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$',
		'datdir3dpp':
			'/home/rpb/compbio/membrane-v612-enthx1-12800/a4-stress-t4-lonestar-75000-175000-200/results',
		'nprots':1,
		'framewise_part':7,
		'voxelsize':1.0,
		'custom_topogcorr_specs':None,
		'whichframes':slice(None,None)},	
	'v701-60000-160000-200':
		{'sysname':'membrane-v701',
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc',
		'nprots':2,
		'datdir3dpp':
			'/home/rpb/compbio-alt/membrane-v701-exo70-anti-dilute/a1-stress-1.0-framewise/results',
		'framewise_part':3,
		'voxelsize':1.0,
		'label':r'$\textbf{{EXO70}\ensuremath{\times}2{\small (anti)}}$',
		'structure_pkl':'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl'},
	'v550-300000-400000-200': 
		{'sysname':'membrane-v550',
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'nprots':0,
		'datdir3dpp':
			'/home/rpb/compbio-alt/membrane-v550/'+\
			'a1-stress-1.0-framewise-md.part0006.300000-400000-200/results',
		'framewise_part':6,
		'voxelsize':1.0,
		'label':r'$\mathrm{control}$'}}
analysis_names = [
	'v614-120000-220000-200',
	'v612-75000-175000-200',
	'v550-300000-400000-200',
	'v614-40000-140000-200'
	][:-1]

#---methods
routine = [
	'calc_c0maps',
	'plot',
	'video',
	'plot_extendview',
	'plot_histograms',
	'plot_2dhistograms'][-2:-1]
span_sweep = [1,2,3,4,5,6]
flatz0 = True #---use a flat plane z0=0 instead of z0=z(x,y) for the integral
bigname = '.'.join(analysis_names)
do_blur = False #---use Gaussian blur
do_blur_first = False #---use Gaussian blur before taking the mean (has no effect)
smooth_wid = 2 #---sigma for the Gaussian blur
window_half = 140 #---half-width in Angstroms of the zoom window for the plot_extendview

#---units
#---Nb this assumes pressure units are in bar, and hence this leads to a factor of ten difference
kappa = 20
temperature = 310
conv_fac = 1.38*10.*kappa*temperature

#---methods
plot_menu = [
	'hist_sum_norm',
	'hist_smooth_sum_norm',
	'join_systems_histogram',
	'systems_join_span_histogram',
	'systems_join_span_histogram_no_mean',
	'neighborhood_unfiltered'
	][-1:]
smooth_mean = 1 #---method: histogram > smooth > sum frames > norm (r)
plot_span = 2 #---method: join all C0 for one span, all systems > histogram
systems_join_span_histogram_subset = False #---method: join all C0 for one span, all systems > histogram

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def stressmap_panel_plot(dat,fig,fr=None,cmap=None,vmax=None,vmin=None,altdat=None,extrarow=False):
	'''Function which plots a stressmap with proteins.'''
	#---settings
	panels = len(dat)
	if cmap == None: cmap = mpl.cm.RdBu_r
	if fr == None: fr = 0
	#---axes
	axeslist = []
	gs = gridspec.GridSpec((2 if extrarow else 1),panels,wspace=0.0,hspace=0.0)
	for p in range(panels):
		print p
		if extrarow: ax = fig.add_subplot(gs[1,p])
		else: ax = fig.add_subplot(gs[p])
		axeslist.append(ax)
		im = ax.imshow((dat[p]).T,interpolation='nearest',origin='lower',
			cmap=cmap,vmax=vmax,vmin=vmin)
		if altdat != None and altdat[p].protein != []:
			plothull(ax,msets[p].protein[fr],griddims=shape(md_maps[0][0].data)[1:],
				vecs=mean(msets[p].vecs,axis=0),subdivide=nnprots[p],alpha=0.35,c='k')
		ax.set_xlabel(r'$x\:(\mathrm{nm})$',fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		#---added the following for flush plots
		if p > 0: ax.set_yticklabels([])
		if p == 0: ax.set_ylabel(r'$y\:(\mathrm{nm})$',fontsize=fsaxlabel)
		if 0: ax.set_title(labels[p],fontsize=fsaxtitle)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
	plt.setp(axins.get_xticklabels(),fontsize=fsaxlabel)
	axins.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel,rotation=270)
	return [gs,axeslist]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'calc_c0maps' in routine or ('md_maps' not in globals() and 'calc_c0maps' not in routine):
	md_maps = [[] for aname in analysis_names]
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
		picklename = 'pkl.stressmaps.'+specname_pickle(sysname,trajfile[0])+\
			('.flat' if flatz0 else '')+'.pkl'
		md_map = unpickle(pickles+picklename)
		if md_map != None: md_maps[analysis_names.index(aname)] = md_map
		labels.append(label)
		nnprots.append(nprots)
		if md_maps[analysis_names.index(aname)] == []:
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
				#---Nb append a suffix to the pickles and use the following line to do the flat calculation
				#ref_surf = mset.unzipgrid(mset.surf[frame],vecs=vecs)
				ref_surf = mset.unzipgrid(zeros(griddims[:2]),vecs=vecs)
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
		
#---a single panel plot
if 'plot' in routine:
	#---average C0 map with PBC Gaussian blur
	m,n = shape(md_maps[0][0].data)[1:]
	#---the following selects the span
	panelplots = [scipy.ndimage.filters.gaussian_filter(numpy.tile(mean(md_maps[i][plot_span].data,axis=0),
		(3,3)),smooth_wid)[m:2*m,n:2*n] for i in range(len(md_maps))]
	vmin = array(panelplots).min()
	vmax = array(panelplots).max()
	extrem = max(abs(vmax),abs(vmin))
	vmax = extrem
	vmin = -1*extrem
	fig = plt.figure()
	gs = stressmap_panel_plot(panelplots,fig,vmax=vmax,vmin=vmin,altdat=msets)
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.6,w_pad=0.6)
	plt.savefig(pickles+'fig-stressmap-'+bigname+'.png',dpi=300,bbox_inches='tight')
	plt.show()
	
#---plot with zoom
if 'plot_extendview' in routine:
	#---average C0 map with PBC Gaussian blur
	m,n = shape(md_maps[0][0].data)[1:]
	#---Nb the Gaussian blur step means that the peaks on the stress maps will be duller than the extrema 
	#---...of the 1D histograms of the average map given in the systems_join_span_histogram step
	if do_blur:
		panelplots = [
			scipy.ndimage.filters.gaussian_filter(numpy.tile(
				mean(md_maps[i][plot_span].data,axis=0),(3,3)),smooth_wid)[m:2*m,n:2*n] 
				for i in range(len(md_maps))]
	elif do_blur_first:
		panelplots = [mean([scipy.ndimage.filters.gaussian_filter(numpy.tile(
			md_maps[i][plot_span].data[j],(3,3)),smooth_wid)[m:2*m,n:2*n] 
			for j in range(len(md_maps[i][plot_span].data))],axis=0) for i in range(len(md_maps))]
	else:
		panelplots = [mean(md_maps[i][plot_span].data,axis=0) for i in range(len(md_maps))]
	#---scale to the correct units
	panelplots = [array(i)/conv_fac for i in panelplots]
	vmin = array(panelplots).min()
	vmax = array(panelplots).max()
	extrem = max(abs(vmax),abs(vmin))
	vmax = extrem
	vmin = -1*extrem
	fig = plt.figure()
	gs,axeslist = stressmap_panel_plot(panelplots,fig,vmax=vmax,vmin=vmin,altdat=msets,extrarow=True)
	#---now draw the extra row with zoomed views on top
	#---Nb protcenter is calculated from msets[0]
	#---get zoomed-in limits here
	protcenter = mean(mean(msets[0].protein,axis=0),axis=0)[:2]
	lims = array([protcenter-array([window_half,window_half]),
		protcenter+array([window_half,window_half])]).T/10.
	vmin = array(panelplots)[:,slice(int(round(lims[0][0])),
		int(round(lims[0][1]))),slice(int(round(lims[0][0])),int(round(lims[0][1])))].min()
	vmax = array(panelplots)[:,slice(int(round(lims[0][0])),
		int(round(lims[0][1]))),slice(int(round(lims[0][0])),int(round(lims[0][1])))].max()
	extrem = max(abs(vmax),abs(vmin))
	vmax = extrem
	vmin = -1*extrem
	#---basically copied the stressmap_panel_plot code right here
	dat = panelplots
	cmap = None
	fr = None
	altdat=msets
	panels = len(dat)
	if cmap == None: cmap = mpl.cm.RdBu_r
	if fr == None: fr = 0
	for p in range(panels):
		print p
		ax = fig.add_subplot(gs[0,p])
		im = ax.imshow((dat[p]).T,interpolation='nearest',origin='lower',
			cmap=cmap,vmax=vmax,vmin=vmin)
		if altdat != None and altdat[p].protein != []:
			plothull(ax,msets[p].protein[fr],griddims=shape(md_maps[0][0].data)[1:],
				vecs=mean(msets[p].vecs,axis=0),subdivide=nnprots[p],alpha=0.2,c='k')
		ax.set_title(labels[p],fontsize=fsaxtitle)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		mset = msets[p] if msets[p].protein != [] else msets[p-1]
		protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
		lims = array([protcenter-array([window_half,window_half]),
			protcenter+array([window_half,window_half])]).T/10.
		ax.set_xlim(lims[0])
		ax.set_ylim(lims[1])
		#---draw a box on the zoomed-out row
		axeslist[p].axvline(x=lims[0][0]-1,ymin=lims[1][0]/(n-0),ymax=lims[1][1]/(n+0),linewidth=2,c='k')
		axeslist[p].axvline(x=lims[0][1],ymin=lims[1][0]/(n-0),ymax=lims[1][1]/(n+0),linewidth=2,c='k')
		axeslist[p].axhline(y=lims[1][0]-1,xmin=lims[0][0]/(m-0),xmax=lims[0][1]/(m+0),linewidth=2,c='k')
		axeslist[p].axhline(y=lims[1][1],xmin=lims[0][0]/(m-0),xmax=lims[0][1]/(m+0),linewidth=2,c='k')
		#---remove for flush plots
		if p > 0: ax.set_yticklabels([])
		if p == 0: ax.set_ylabel(r'$y\:(\mathrm{nm})$',fontsize=fsaxlabel)
		if 0: ax.set_xlabel(r'$x\:(\mathrm{nm})$',fontsize=fsaxlabel)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
	plt.setp(axins.get_xticklabels(),fontsize=fsaxlabel)
	axins.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel,rotation=270)
	#---save
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.0,w_pad=0.0)
	plt.savefig(pickles+'fig-stressmap_zoom-'+bigname+('.flat' if flatz0 else '')+\
		'-span'+str(span_sweep[plot_span])+\
		('-blur'+str(smooth_wid) if do_blur else '')+\
		('-blurfirst'+str(smooth_wid) if do_blur_first else '')+\
		'.png',dpi=300,bbox_inches='tight')
	plt.show()
	
#---make cool videos
if 'video' in routine:
	vidpanels = []
	for p in range(len(md_maps)):
		vidpanels.append(array([scipy.ndimage.filters.gaussian_filter(mean(md_maps[p][plot_span].data[i:i+20],
			axis=0),2) for i in range(500-20)]))
	plotmov(vidpanels,'stressmap-'+bigname,altdat=msets,panels=len(md_maps),
		plotfunc='stressmap_panel_plot',whitezero=True)

#---extra list of possible histograms
for plot_type in plot_menu:
	#---method: histogram > sum frames > norm (r)	
	if plot_type == 'hist_sum_norm':
		fig = plt.figure()
		ax = plt.subplot(111)
		#---sum the histogram over protein distance and C0 magnitude then normalize by protein distance
		dat = sum(bigdathist,axis=0).T/sum(sum(bigdathist,axis=0),axis=1)
		ax.imshow(dat,interpolation='nearest',origin='lower')
		plt.show()
	#---method: histogram > smooth > sum frames > norm (r)	
	if plot_type == 'hist_smooth_sum_norm':
		#---compute
		smoothdat = [scipy.ndimage.filters.gaussian_filter(i,smooth_mean) for i in bigdathist]
		summeddat = sum(smoothdat,axis=0).T/sum(sum(smoothdat,axis=0),axis=1)
		#---plot
		fig = plt.figure()
		ax = plt.subplot(111)
		ax.imshow(summeddat,interpolation='nearest',origin='lower',
			cmap=mpl.cm.jet)
		plt.plot([sum(summeddat[:,i]*range(60))/sum(summeddat[:,i]) for i in range(70)],c='k',lw=3)
		plt.show()
	#---method: for each system join c0 for all frames for all spans > histogram
	if plot_type == 'join_systems_histogram':
		fig = plt.figure()
		ax = plt.subplot(111)
		for j in range(len(md_maps)):
			md_map = md_maps[j]
			for i in range(len(md_map)):
				counts,edges = histogram(array(md_map[i].data).flatten())
				mids = (edges[1:]+edges[:-1])/2.
				plt.plot(mids,counts,lw=2,c=clrs[j])
		plt.show()
	#---method: for each system join c0 for all frames for one span > histogram
	if plot_type == 'systems_join_span_histogram':
		width = 0.35
		nbins = 21 if systems_join_span_histogram_subset else 41
		fig = plt.figure()
		ax = plt.subplot(111)
		ax.axvline(x=0,ymin=0,ymax=1,color='k',lw=2)
		pozsums = []
		negsums = []
		for j in range(len(md_maps)):
			md_map = md_maps[j]
			#---box subset
			if systems_join_span_histogram_subset:
				mset = msets[j] if msets[j].protein != [] else msets[j-1]
				protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
				boxlims = array([protcenter-array([window_half,window_half]),
					protcenter+array([window_half,window_half])]).T/10.
				boxlims = [[int(round(i)) for i in l] for l in boxlims]
				dat = mean(array(array(md_map[plot_span].data)[:,boxlims[0][0]:boxlims[0][1],
					boxlims[1][0]:boxlims[1][1]]),axis=0).flatten()/conv_fac
			#---full box
			else: dat = mean(md_map[plot_span].data,axis=0).flatten()/conv_fac
			lims = [[-1*i,i] for i in [max([abs(dat.min()),abs(dat.max())])]][0]
			counts,edges = histogram(dat,range=lims,bins=nbins)
			mids = (edges[1:]+edges[:-1])/2.
			ax.plot(mids,counts,'o-',lw=2,c=clrs[j],
				label=(analysis_descriptors[analysis_names[j]])['label'],mec=clrs[j])
			pozsums.append(sum(counts[mids>=0]*mids[mids>=0]))
			negsums.append(abs(sum(counts[mids<=0]*mids[mids<=0])))
		#---inset comparing positive and negative values
		axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="25%",height="40%",loc=1)
		ind = np.arange(len(md_maps))
		pozbars = axins.bar(ind+width,pozsums,width,color='r',alpha=0.5)
		negbars = axins.bar(ind,negsums,width, color='b',alpha=0.5)
		axins.set_xticklabels([(analysis_descriptors[analysis_names[j]])['label'] 
			for j in range(len(analysis_names))])
		plt.setp(axins.get_xticklabels(), rotation=90)
		axins.set_yticklabels([])
		axins.set_xticks(ind+width)
		axins.set_ylabel(r'${C}_{0}\,\mathrm{balance}$',fontsize=fsaxlabel)
		axins.set_yticks([])
		plt.tick_params(axis='x',which='both',top='off')
		#---plot
		ax.set_xlabel(r'$\mathsf{C_{0}(x,y)\,(nm^{-1})}$',fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		ax.set_yticklabels([])
		ax.grid(True)
		ax.legend(loc='upper left')
		plt.savefig(pickles+'fig-stressdist-joinhist-'+\
			bigname+'-span'+str(span_sweep[plot_span])+\
			('.flat' if flatz0 else '')+\
			('.subset' if systems_join_span_histogram_subset else '')+\
			'.png',dpi=300,bbox_inches='tight')
		plt.show()
	#---method: for each system join c0 for all frames with no mean for one span > histogram
	if plot_type == 'systems_join_span_histogram_no_mean':
		width = 0.35
		nbins = 21 if systems_join_span_histogram_subset else 41
		fig = plt.figure()
		ax = plt.subplot(111)
		ax.axvline(x=0,ymin=0,ymax=1,color='k',lw=2)
		pozsums = []
		negsums = []
		extrem = 0
		for j in range(len(md_maps)):
			md_map = md_maps[j]
			#---box subset
			if systems_join_span_histogram_subset:
				mset = msets[j] if msets[j].protein != [] else msets[j-1]
				protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
				boxlims = array([protcenter-array([window_half,window_half]),
					protcenter+array([window_half,window_half])]).T/10.
				boxlims = [[int(round(i)) for i in l] for l in boxlims]
				dat = array(array(md_map[plot_span].data)[:,boxlims[0][0]:boxlims[0][1],
					boxlims[1][0]:boxlims[1][1]]).flatten()/conv_fac
			#---full box
			else: dat = array(md_map[plot_span].data).flatten()/conv_fac
			thismax = max([abs(dat.min()),dat.max()])
			if thismax > extrem: extrem = thismax
		lims = [-thismax,thismax]
		for j in range(len(md_maps)):
			md_map = md_maps[j]
			#---box subset
			if systems_join_span_histogram_subset:
				mset = msets[j] if msets[j].protein != [] else msets[j-1]
				protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
				boxlims = array([protcenter-array([window_half,window_half]),
					protcenter+array([window_half,window_half])]).T/10.
				boxlims = [[int(round(i)) for i in l] for l in boxlims]
				dat = array(array(md_map[plot_span].data)[:,boxlims[0][0]:boxlims[0][1],
					boxlims[1][0]:boxlims[1][1]]).flatten()/conv_fac
			#---full box
			else: dat = array(md_map[plot_span].data).flatten()/conv_fac
			#lims = [[-1*i,i] for i in [max([abs(dat.min()),abs(dat.max())])]][0]
			counts,edges = histogram(dat,range=lims,bins=nbins)
			mids = (edges[1:]+edges[:-1])/2.
			ax.plot(mids,counts,'o-',lw=2,c=clrs[j],
				label=(analysis_descriptors[analysis_names[j]])['label'],mec=clrs[j])
			pozsums.append(sum(counts[mids>=0]*mids[mids>=0]))
			negsums.append(abs(sum(counts[mids<=0]*mids[mids<=0])))
		#---inset comparing positive and negative values
		axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="25%",height="40%",loc=1)
		ind = np.arange(len(md_maps))
		pozbars = axins.bar(ind+width,pozsums,width,color='r',alpha=0.5)
		negbars = axins.bar(ind,negsums,width, color='b',alpha=0.5)
		axins.set_xticklabels([(analysis_descriptors[analysis_names[j]])['label'] 
			for j in range(len(analysis_names))])
		plt.setp(axins.get_xticklabels(), rotation=90)
		axins.set_yticklabels([])
		axins.set_xticks(ind+width)
		axins.set_ylabel(r'${C}_{0}\,\mathrm{balance}$',fontsize=fsaxlabel)
		axins.set_yticks([])
		plt.tick_params(axis='x',which='both',top='off')
		#---plot
		ax.set_xlabel(r'$\mathsf{C_{0}(x,y)\,(nm^{-1})}$',fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		ax.set_yticklabels([])
		ax.grid(True)
		ax.legend(loc='upper left')
		plt.savefig(pickles+'fig-stressdist-joinhistall-'+\
			bigname+'-span'+str(span_sweep[plot_span])+\
			('.flat' if flatz0 else '')+\
			('.subset' if systems_join_span_histogram_subset else '')+\
			'.png',dpi=300,bbox_inches='tight')
		plt.show()
	#---method: raw histograms in three different domains
	if plot_type == 'neighborhood_unfiltered':
		'''
		Notes: in the section called 'systems_join_span_histogram_no_mean', we can plot the raw data in a 
		...histogram for the full box and subset. In this section, I'm rolling this into a panel plot instead 
		...of having to switch the flags and do separate figures. I'm also adding a third method which allows
		...you to do a raw histogram of points within a certain distance of the protein.
		'''
		width = 0.35
		nbins = 101
		fig = plt.figure(figsize=(6,10))
		gs = gridspec.GridSpec(len(msets),1,wspace=0.0,hspace=0.0)
		extrem = 0
		boxslice = ['full','center','near']
		marklist = ['o','v','>']
		#---first find the correct limits
		for j in range(len(md_maps)):
			for k in boxslice:
				md_map = md_maps[j]
				if k == 'center':
					mset = msets[j] if msets[j].protein != [] else msets[j-1]
					protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
					boxlims = array([protcenter-array([window_half,window_half]),
						protcenter+array([window_half,window_half])]).T/10.
					boxlims = [[int(round(i)) for i in l] for l in boxlims]
					dat = array(array(md_map[plot_span].data)[:,boxlims[0][0]:boxlims[0][1],
						boxlims[1][0]:boxlims[1][1]]).flatten()/conv_fac
				elif k == 'full': dat = array(md_map[plot_span].data).flatten()/conv_fac
				elif k == 'near':
					#---take the mean protein position and identify all points on the md_maps grid that are 
					#---...within some cutoff of the time-wise averaged protein locations then use these 
					#---...points in the histogram
					m,n = msets[j].griddims
					protpts = mean(msets[0].protein,axis=0)
					mvecs = mean(msets[0].vecs,axis=0)
					mapper = md_maps[j][plot_span].data[0]
					surfpts = mset.wrappbc(mset.unzipgrid(array(mapper),vecs=mset.vec(0)),vecs=mvecs)
					cd = scipy.spatial.distance.cdist(protpts[:,:2],surfpts[:,:2])
					dists = array([np.min(cd,axis=0),surfpts[:,2]]).T
					shadow_inds = where(np.min(cd,axis=0)<=200)[0]
					shadow_inds = array(unravel_index(shadow_inds,shape(mapper))).T
					dat = array([[md_map[plot_span].data[fr][i[0],i[1]] for i in shadow_inds] 
						for fr in range(len(md_map[plot_span].data))]).flatten()/conv_fac
				thismax = max([abs(dat.min()),dat.max()])
				if thismax > extrem: extrem = thismax
		#---perform the plotting
		lims = [-thismax,thismax]
		for j in range(len(md_maps)):
			ax = fig.add_subplot(gs[j])
			ax.axvline(x=0,ymin=0,ymax=1,color='k',lw=2)
			pozsums = []
			negsums = []
			for k in boxslice:
				md_map = md_maps[j]
				if k == 'center':
					mset = msets[j] if msets[j].protein != [] else msets[j-1]
					protcenter = mean(mean(mset.protein,axis=0),axis=0)[:2]
					boxlims = array([protcenter-array([window_half,window_half]),
						protcenter+array([window_half,window_half])]).T/10.
					boxlims = [[int(round(i)) for i in l] for l in boxlims]
					dat = array(array(md_map[plot_span].data)[:,boxlims[0][0]:boxlims[0][1],
						boxlims[1][0]:boxlims[1][1]]).flatten()/conv_fac
				elif k == 'full': dat = array(md_map[plot_span].data).flatten()/conv_fac
				elif k == 'within 10 nm':
					#---take the mean protein position and identify all points on the md_maps grid that are 
					#---...within some cutoff of the time-wise averaged protein locations then use these 
					#---...points in the histogram
					m,n = msets[j].griddims
					protpts = mean(msets[0].protein,axis=0)
					mvecs = mean(msets[0].vecs,axis=0)
					mapper = md_maps[j][plot_span].data[0]
					surfpts = mset.wrappbc(mset.unzipgrid(array(mapper),vecs=mset.vec(0)),vecs=mvecs)
					cd = scipy.spatial.distance.cdist(protpts[:,:2],surfpts[:,:2])
					dists = array([np.min(cd,axis=0),surfpts[:,2]]).T
					shadow_inds = where(np.min(cd,axis=0)<=200)[0]
					shadow_inds = array(unravel_index(shadow_inds,shape(mapper))).T
					dat = array([[md_map[plot_span].data[fr][i[0],i[1]] for i in shadow_inds] for 
						fr in range(len(md_map[plot_span].data))]).flatten()/conv_fac
				
				counts,edges = histogram(dat,range=lims,bins=nbins,normed=True)
				mids = (edges[1:]+edges[:-1])/2.
				if 0: ax.plot(mids,counts,'o-',lw=1,c=clrs[boxslice.index(k)],
					label=k,mec=clrs[boxslice.index(k)],
					marker=marklist[boxslice.index(k)])

				if 1: ax.plot(mids,counts,'o-',lw=1,c=clrs[boxslice.index(k)],
					label=k)


				pozsums.append(sum(counts[mids>=0]*mids[mids>=0]))
				negsums.append(abs(sum(counts[mids<=0]*mids[mids<=0])))
				
			if 1:
				ax.set_xlim((-0.05,0.05))
				#ax.set_ylim((5.,6.0))
				ax.set_ylim((0.8*max(counts),max(counts)*1.1))
			
			ax.set_ylabel((analysis_descriptors[analysis_names[j]])['label'],fontsize=fsaxlabel)
			if j < len(msets)-1: ax.set_xticklabels([])

			#---inset comparing positive and negative values
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="25%",height="40%",loc=1)
			ind = np.arange(len(md_maps))
			pozbars = axins.bar(ind+width,pozsums,width,color='r',alpha=0.5)
			negbars = axins.bar(ind,negsums,width, color='b',alpha=0.5)
			axins.set_xticklabels([(analysis_descriptors[analysis_names[j]])['label'] 
				for j in range(len(analysis_names))])
			plt.setp(axins.get_xticklabels(), rotation=90)
			axins.set_yticklabels([])
			axins.set_xticks(ind+width)
			axins.set_ylabel(r'${C}_{0}\,\mathrm{balance}$',fontsize=fsaxlabel)
			axins.set_yticks([])
			plt.tick_params(axis='x',which='both',top='off')
			#---plot
			ax.set_xlabel(r'$\mathsf{C_{0}(x,y)\,(nm^{-1})}$',fontsize=fsaxlabel)
			plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
			ax.set_yticklabels([])
			ax.grid(True)
			ax.legend(loc='upper left')
		'''
		plt.savefig(pickles+'fig-stressdist-joinhistall-'+\
			bigname+'-span'+str(span_sweep[plot_span])+\
			('.flat' if flatz0 else '')+\
			('.subset' if systems_join_span_histogram_subset else '')+\
			'.png',dpi=300,bbox_inches='tight')
		'''
		plt.show()
	
#---plot 2D histograms over protein distance and and curvature
if 'plot_2dhistograms' in routine and 'stack_collected_hists' not in globals():
	nbins = 20
	sourcedatflat = array(md_maps[0][plot_span].data).flatten()
	nhistmax = max([abs(sourcedatflat.min()),abs(sourcedatflat.max())])
	binw = (10,100)
	stack_collected_dist_vs_c0 = []
	stack_collected_hists = []
	for j in range(len(md_maps)):
		mset = msets[j]
		mapper = md_maps[j][plot_span].data
		protsref = slice(None,None)
		cutoff = min(mean(mset.vecs,axis=0)[:2])
		hist_ranges = ((0,ceil(cutoff/100)*100),(-nhistmax,nhistmax))
		nbins = [round((hist_ranges[i][1]-hist_ranges[i][0])/binw[i]) for i in range(2)]
		nbins = [nbins[0],(nbins[1]+((nbins[1]-1)%2))]
		bigdat = []
		bigdathist = []
		#---if no protein use the previous panel
		mset_protein = msets[j] if msets[j].protein != [] else msets[j-1]
		for fr in range(len(mapper)):
			print analysis_names[j]+' '+str(fr)
			protpts = mset_protein.protein[fr][protsref,0:2]
			surfpts = mset.wrappbc(mset.unzipgrid(array(mapper[fr]),vecs=mset.vec(0)),
				vecs=mset.vec(fr),mode='nine')
			#---Nb torus_metric works here
			#---Nb that the torus metric works fine for distance
			#---Nb but the torus metric doesn't replicate points the way we want
			cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2])
			bigdat.append(array([np.min(cd,axis=0),surfpts[:,2]]).T)
			H, xedges, yedges = histogram2d(bigdat[-1][:,0],bigdat[-1][:,1],bins=nbins,range=hist_ranges)
			#row_sums = H.sum(axis=1)
			bigdathist.append(array(H))
		stack_collected_dist_vs_c0.append(bigdat)
		stack_collected_hists.append(bigdathist)
if 'plot_2dhistograms' in routine:
	aspect = 1./3.
	binspan = 10
	fig = plt.figure(figsize=(6,10))
	gs = gridspec.GridSpec(1,len(md_maps),wspace=0.0,hspace=0.0)
	if 0: ax = fig.add_subplot(gs[p])
	if 0: ax = plt.subplot(111)
	for p in range(len(md_maps)):
		ax = fig.add_subplot(gs[p])
		sourcedat = stack_collected_hists[p]
		m,n = shape(sourcedat)[1:]
		dat = mean(sourcedat,axis=0).T/mean(mean(sourcedat,axis=0),axis=1)
		im = ax.imshow(dat,interpolation='nearest',origin='lower',aspect=2.,
			vmax=dat.max()*1.1,vmin=dat.min(),cmap=mpl.cm.binary,alpha=0.5)
		if 0: extent =  im.get_extent()
		if 0: ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
		means = array([sum(dat[:,i]*range(n))/sum(dat[:,i]) for i in range(m)])
		stddevs = array([sqrt(sum(dat[:,i]/sum(dat[:,i])*(range(len(dat))-means[i])**2)) for i in range(m)])
		ax.axhline(y=int(nbins[1]/2),xmin=0.,xmax=1.,color='w',lw=3)
		ax.plot([sum(dat[:,i]*range(n))/sum(dat[:,i]) for i in range(m)],lw=4,color='k',
			label=(analysis_descriptors[analysis_names[p]])['label'])
		if 0: ax.fill_between(range(m),means+stddevs,means-stddevs,color='k',alpha=0.2)
		ax.set_ylim((nbins[1]/2.-binspan,nbins[1]/2.+binspan))
		yticks = [int(round(i)) for i in linspace(n/2-binspan,n/2+binspan,7)]
		ylabels = [round(int(round(linspace(-nhistmax,nhistmax,n)[j],-2))/conv_fac,3) for j in yticks]
		ax.set_yticklabels(ylabels)
		ax.set_yticks(yticks)
		ax.set_xlim((0,20))
		if p > 0:
			ax.set_yticklabels([])
		if p == 0: ax.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel)
		ax.set_title((analysis_descriptors[analysis_names[p]])['label'],fontsize=fsaxlabel)
		ax.set_xlabel(r'$\left|\mathbf{r}_{min}\right|\,\mathrm{(nm)}$',fontsize=fsaxlabel)
		ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper',nbins=5))
	if 0: plt.legend()
	gs.tight_layout(fig,h_pad=0.0,w_pad=0.0)
	plt.savefig(pickles+'fig-stressdist-distvc0-'+\
		bigname+'-span'+str(span_sweep[plot_span])+\
		('.flat' if flatz0 else '')+\
		'.png',dpi=300,bbox_inches='tight')
	plt.show()
	

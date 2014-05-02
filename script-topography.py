#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

import scipy.ndimage
from scipy.signal import hilbert
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt
from scipy.spatial.distance import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid import Size
from mpl_toolkits.axes_grid import Divider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---plan
analysis_descriptors_extra = {
	'v614-120000-220000-200': {
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'topogcorr_pkl':'pkl.topogcorr.membrane-v614-s9-lonestar-120000-220000-200.pkl'
		},
	'v612-75000-175000-200': {
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'whichframes':slice(None,None)
		},
	'v550-400000-500000-160': {
		'protein_pkl':'pkl.structures.space10A.membrane-v612.t4-lonestar.md.part0007.75000-175000-200.pkl',
		'custom_topogcorr_specs':None,
		'custom_protein_shifts':['unshifted','peak','valley',[0.,-70.]][1:-1],
		},
	}
		
#---access central dictionary and combine it with the specifications above
execfile('header-cgmd.py')
		
#---analysis menu
do = ['phasecomp','topography_calc','topography_plot'][2]
if do == 'phasecomp':
	analysis_names = [
		'v700-500000-600000-200',
		'v614-120000-220000-200',
		'v550-400000-500000-160'
		][:]
	plot_reord = analysis_names
	dolist = ['phase_std_spec2d']
	bigname = '.'.join(analysis_names)
	nhistbins = 80
elif do == 'topography_calc':
	analysis_names = [
		'v614-120000-220000-200',
		'v612-75000-175000-200',
		'v550-400000-500000-160'
		][-1:]
	plot_reord = analysis_names
	dolist = ['topogcorr_calc']
	nhistbins = 80
elif do == 'topography_plot':
	analysis_names = [
		'v614-120000-220000-200',
		'v612-75000-175000-200',
		'v550-400000-500000-160'
		][:]
	plot_reord = analysis_names
	dolist = ['topogcorr_calc','topogcorr_plot_alt'][1:]
	topogcorr_plot_zrange = (-2.0,2.0)
	bigname = 'v614-v612-v550-ver2'
	nhistbins = 80
	nhistbins_smooth = 10

#---settings
cmap = mpl.cm.RdBu_r
clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in range(9)]
clrdict = {'red':0,'blue':1}

#---method, topographic correlate
topogcorr_plot_struct = True

#---specify grid spacing for more recent midplane data
gridspacing = 1.0
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'

#---specify how many frames to use
whichframes = slice(0,500)

#---MAIN, PHASE PLOTS
#-------------------------------------------------------------------------------------------------------------

def specialname(basename):
	'''Standard method for deriving a step/part/timeslice-specific name from the trajectory file.'''
	special_name = "-".join([i for i in re.match('.+membrane\-.+',basename).string.split('.')
		if re.match('[a-z][0-9]\-.+',i) 
		or re.match('membrane-v.+',i) 
		or re.match('[0-9]+\-[0-9]+\-[0-9]',i)])
	return special_name
	
def torus_metric(x1,x2,vecs,ndims=2):
	'''The torus metric i.e. how to measure distance in periodic spaces, quickly.'''
	raw = array([cdist(x1[:,i:i+1],x2[:,i:i+1]) for i in range(ndims)])
	pbcd = array([raw[i]-(raw[i]>vecs[i]/2)*vecs[i] for i in range(ndims)]).T
	return sqrt(sum(pbcd**2,axis=-1)).T

#---MAIN, PHASE PLOTS
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		msets.append(unpickle(pickles+'pkl.structures.'+spacetag+specname_guess(sysname,trajsel)+'.pkl'))
		msets[-1].calculate_undulations()
		#---hack
		msets[-1].rounder = 20.

#---phase plots
#---Nb this is a canonical spectrum-plotting function, consider standardizing and adding to plotter.py
if 'phase_std_spec2d' in dolist:
	calcs = ['mean','var']
	calcnamesl = ['mean phase','phase variance']
	calcnamesr = [r'$\mathrm{rad}$',r'$\mathrm{rad}^{2}$']
	smooth_std = 1
	smooth_mean = 1
	#---prepare data	
	results_dat = [[[] for i in range(len(analysis_names))] for j in range(len(calcs))]
	for calc in calcs:
		for m in [analysis_names.index(aname) for aname	in plot_reord]:
			a = analysis_names[m]
			for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
			mset = msets[m]
			if calc == 'var':
				allangles = array([angle(mset.undulate_raw[i]) for i in range(len(mset.undulate_raw))])
				if smooth_std != None:
					results_dat[calcs.index(calc)][m] = \
						scipy.ndimage.filters.gaussian_filter(var(allangles,axis=0),smooth_std)
				else:
					results_dat[calcs.index(calc)][m] = var(allangles,axis=0)
			elif calc == 'mean':
				allangles = array([angle(mset.undulate_raw[i]) for i in range(len(mset.undulate_raw))])
				if smooth_mean != None:
					results_dat[calcs.index(calc)][m] = \
						scipy.ndimage.filters.gaussian_filter(mean(allangles,axis=0),smooth_mean)
				else:
					results_dat[calcs.index(calc)][m] = mean(allangles,axis=0)
	#---prepare figure
	fig = plt.figure()
	#---Nb using the axesgrid objects here, despite the restriction that they must shrae certain axes
	#---Nb the axesgrid object lets you make flush equal-sized axes with flush colorbars
	gs = [AxesGrid(fig,211,nrows_ncols=(1,len(analysis_names)),axes_pad=0.0,share_all=False,label_mode='L',
		cbar_location="right",cbar_mode="single"),
		AxesGrid(fig,212,nrows_ncols=(1,len(analysis_names)),axes_pad=0.0,share_all=False,label_mode='L',
        cbar_location="right",cbar_mode="single")]
	for calc in calcs:
		c = calcs.index(calc)
		vmax = array(results_dat[c]).max()
		vmin = array(results_dat[c]).min()
		#---plot undulations
		for m in [analysis_names.index(aname) for aname	in plot_reord]:
			a = analysis_names[m]
			for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
			mset = msets[m]
			data = results_dat[calcs.index(calc)][m]
			#---plot
			#ax = plt.subplot(gs[calcs.index(calc),(plot_reord).index(a)])
			ax = gs[calcs.index(calc)][(plot_reord).index(a)]
			im = plotter2d(ax,mset,dat=data,cmap=cmap,lognorm=False,lims=[vmin,vmax],
				ticklabel_show=[(1 if c == len(calcs)-1 else 0),0],tickshow=True,label_style='q',
				centertick=True,tickskip=int(round(mset.griddims[0]/6,-1)))
			#---plot details
			if c == 0: ax.set_title(label)
		#---colorbar for this row
		cbar = gs[c].cbar_axes[0].colorbar(im)
		gs[c].cbar_axes[0].set_ylabel(calcnamesr[c],fontsize=fsaxlabel)
		gs[c][0].set_ylabel(calcnamesl[c],fontsize=fsaxlabel)
	plt.savefig(pickles+'fig-phases-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.show()

#---MAIN, TOPOGRAPHY
#-------------------------------------------------------------------------------------------------------------

class TopographyCorrelate():
	'''Class object for holding analysis of topography correlation data.'''
	def __init__(self,mset,nprots,label,mset_protein=None):
		self.histdat = []
		self.nprots = nprots
		#---optional use of another source of protein points
		if mset_protein == None: mset_protein = mset
		self.prots_nres = shape(mset_protein.protein[0])[0]/nprots
		#---this class is largely just a wrapper for a list of topographic matrices
		if protein_pkl != None:
			self.prots_refs = []
		else:
			self.prots_refs = ([slice(None,None)] if nprots > 1 else [])+\
				[slice(i*self.prots_nres,(i+1)*self.prots_nres) 
				for i in range(nprots)]
		self.prots_names = [label]+['domain '+str(i) for i in range(nprots)]
		self.prots_shifts = []
		self.cutoff = min(mean(mset.vecs,axis=0)[:2])
		#---list of matrices of protein-to-surface point distances
		self.topogmat = []
		#---list of histogram matrices of protein-to-surface vs surface height distances
		self.topoghist = []
		#---default bin widths
		self.binw = (10,1)
		#---post-processed mean correlation maps
		self.hcorr = []
		self.notes = []
	def makevid(self,index=0):
		'''Makes a video of a single topography histogram, specified by index.'''
		plotmov(tc.topoghist,'test-heightcorr',smooth=1,
			xedges=[int(i) for i in xedges],yedges=[int(i) for i in yedges])
	def addnote(self,notable):
		self.notes.append(notable)

#---topography correlate analysis
if 'topogcorr_calc' in dolist:
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None:
			mset_protein = unpickle(pickles+protein_pkl)
		else:
			mset_protein = mset
		tc = TopographyCorrelate(mset,nprots,label,mset_protein=mset_protein)
		for i in analysis_descriptors[a]: 
			tc.addnote([i,str((analysis_descriptors[a])[i])])
		#---additional protein sub-domains to test
		if custom_topogcorr_specs != None:
			for item in custom_topogcorr_specs:
				tc.prots_refs.append(item[0])
				tc.prots_names.append(item[1])
		#---if studying a control, get protein points from another pickle and shift them
		elif protein_pkl != None:
			for item in custom_protein_shifts:
				#---default is to use the entire set of protein points from the alternate pickle
				tc.prots_refs.append(slice(None,None))
				tc.prots_names.append(item)
				#---we only align the average COM of the protein to the average peak or valley position
				if item == 'peak':
					maxpos = unravel_index(mean(mset.surf,axis=0).argmax(),mset.surf[0].shape)
					maxxy = [maxpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
					protein_com = mean(mean(mset_protein.protein,axis=0),axis=0)[:2]
					shift = maxxy - protein_com
					tc.prots_shifts.append(shift)
				elif item == 'valley':
					maxpos = unravel_index(mean(mset.surf,axis=0).argmin(),mset.surf[0].shape)
					minxy = [maxpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
					protein_com = mean(mean(mset_protein.protein,axis=0),axis=0)[:2]
					shift = minxy - protein_com
					tc.prots_shifts.append(shift)
				elif item == 'unshifted':
					tc.prots_shifts.append(array([0.,0.]))
				elif type(item) == list:
					tc.prots_shifts.append(array(item))
		for pnum in range(len(tc.prots_refs)):
			protsref = tc.prots_refs[pnum]
			hist_ranges = ((0,ceil(tc.cutoff/100)*100),(-nhistbins/2,nhistbins/2))
			nbins = [(hist_ranges[i][1]-hist_ranges[i][0])/tc.binw[i] for i in range(2)]
			topogmat = []
			topoghist = []
			for fr in range(len(mset.surf[whichframes])):
				print 'subset: '+str(protsref)+', frame: '+str(fr)
				if tc.prots_shifts != []:
					shift = tc.prots_shifts[pnum]
					protpts = mset_protein.protein[fr][protsref,0:2]+shift
				else:
					protpts = mset_protein.protein[fr][protsref,0:2]
				surfpts = mset.wrappbc(mset.unzipgrid(array(mset.surf[fr]),vecs=mset.vec(0)),
					vecs=mset.vec(fr),mode='nine')
				#---Nb torus_metric works here
				#---Nb that the torus metric works fine for distance
				#---Nb but the torus metric doesn't replicate points the way we want
				cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2])
				topogmat.append(array([np.min(cd,axis=0),surfpts[:,2]]).T)
				H, xedges, yedges = histogram2d(topogmat[-1][:,0],topogmat[-1][:,1],
					bins=nbins,range=hist_ranges)
				row_sums = H.sum(axis=1)
				topoghist.append(array(H))
			#---Normalize, then take the mean, resulting in a much smoother curve
			histdatnorm = nan_to_num([array(topoghist[i]/(topoghist[i].sum(axis=0))) 
				for i in range(len(topoghist))])
			hcorr = mean(histdatnorm,axis=0)
			tc.topoghist.append(array(topoghist))
			tc.topogmat.append(array(topoghist))
			tc.hcorr.append(array(hcorr))
		tc.xedges = xedges
		tc.yedges = yedges
		#---wasteful to save the entire distance matrix, so recompute later if necessary
		del tc.topogmat
		pickledump(tc,'pkl.topography_correlate.'+specname_guess(sysname,trajsel)+'.pkl',
			directory=pickles)

#---plot the topography correlate curves deprecated by the alternate method below
if 'topogcorr_plot' in dolist:
	zrange = topogcorr_plot_zrange if len(analysis_names) > 1 else None
	swaprows = True if len(analysis_names) == 1 else False
	fig = plt.figure()
	gs = gridspec.GridSpec((len(analysis_names) if swaprows else 2),
		(2 if swaprows else len(analysis_names)),wspace=0.0,hspace=0.0)
	tcs = []
	#---loop over analyses
	for a in analysis_names:
		#---prepare figure with gridspec
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None:
			mset_protein = unpickle(pickles+protein_pkl)
		else:
			mset_protein = mset
		tc = unpickle(pickles+'pkl.topography_correlate.'+specname_guess(sysname,trajsel)+'.pkl')
		tcs.append(tc)
		smoothwid = 4
		ax = fig.add_subplot((gs[0] if len(analysis_names) == 1 
			else gs[(m if swaprows else 0),(0 if swaprows else m)]))
		hcorr_smooths = []
		for h in range(len(tc.hcorr)):
			#---set labels
			if h == 0: color = 'k'
			elif h < len(tc.prots_refs): color = clrs[(h-1)%(len(tc.prots_refs)-1)]
			#---Nb since we use the argmax to give the mode of a broad distribution, hard to find error bars
			#---Nb tried taking the second moment relative to the mode but this made no sense
			histdatnorm = nan_to_num([array(tc.topoghist[h][i]/(tc.topoghist[h][i].sum(axis=0))) 
				for i in range(len(tc.topoghist[h]))])
			hcorr = mean(histdatnorm,axis=0)
			hcorr_smooth = [mean([tc.yedges[hcorr[j].argmax()] for j in range(i+smoothwid)]) 
				for i in range(len(hcorr)-smoothwid)]
			hcorr_smooth_std = [std([tc.yedges[hcorr[j].argmax()] for j in range(i+smoothwid)]) 
				for i in range(len(hcorr)-smoothwid)]
			#---extra attempt to distinguish control from protein simulation			
			if 0:
				extra_smoothwid = 20
				tmp = [10*std(mean(mean(tc.topoghist[0][:,s:s+extra_smoothwid],axis=0),axis=0)) 
					for s in range(0,len(hcorr)-extra_smoothwid,1)]
				hcorr_smooth_err = tmp+(len(hcorr)-len(tmp)-smoothwid)*tmp[-1:]
			else:
				hcorr_smooth_err = hcorr_smooth_std
			#---override names if you have a separate protein pickle
			if protein_pkl == None:
				if h == 0:
					plotlabel = 'full'
				else:
					plotlabel = tc.prots_names[h]
			else:
				plotlabel = tc.prots_names[h+2]
			if protein_pkl == None or tc.prots_names[h+2] in custom_protein_shifts:
				ax.plot(array(tc.xedges[1:-smoothwid])/10.,array(hcorr_smooth)/10.,'o-',label=plotlabel,
					c=color)
				ax.fill_between(array(tc.xedges[1:-smoothwid])/10.,
					array(hcorr_smooth)/10.+array(hcorr_smooth_err)/10.,
					array(hcorr_smooth)/10.-array(hcorr_smooth_err)/10.,
					color=color,alpha=0.1)
				if len(tc.prots_refs) > 1:
					ax.legend(fontsize=fsaxlegend)
				ax.grid(True)
		ax.set_xlabel(r'$r_{min}\:(\mathrm{nm})$',
			fontsize=fsaxlabel)
		ax.set_ylabel(r'$\left\langle z_{mode}\right\rangle (\mathrm{nm})$',fontsize=fsaxlabel)
		ax.set_title(label,fontsize=fsaxtitle)
		ax.axhline(y=0,xmin=0.,xmax=max(tc.xedges[1:-smoothwid]),lw=2,color='k')
		ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper'))
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		ax.set_ylim(zrange)
		#---clean up
		if m != 0:
			ax.set_yticklabels([])
			ax.set_ylabel('')
		#---bottom row shows average structures
		ax = fig.add_subplot((gs[1] if len(analysis_names) == 1 
			else gs[(m if swaprows else 1),(1 if swaprows else m)]))
		maxval = max(abs(mean(mset.surf,axis=0).min()/10.),abs(mean(mset.surf,axis=0).max()/10.))
		im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=cmap,
			inset=False,cmap_washout=1.0,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=fsaxlabel,label_style='xy',
				lims=([-maxval,maxval] if zrange == None else [zrange[0],zrange[1]]))
		s0 = 1 if protein_pkl == None else 0
		s0 = 0 if len(tc.prots_refs) == 1 else s0
		for s in range(s0,len(tc.prots_refs)):
			if protein_pkl == None or tc.prots_names[s+2] in custom_protein_shifts:
				color = 'k' if s == 0 else clrs[s-s0]
				if not hasattr(tc,'prots_shifts'): shift = array([0.,0.])
				elif tc.prots_shifts == []: shift = array([0.,0.])
				else: shift = tc.prots_shifts[s]
				for pbcshift in [[0,0],[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]:
					plothull(ax,mean(mset_protein.protein,axis=0)[tc.prots_refs[s],:2]+\
						shift+array([mean(mset.vecs,axis=0)[i]*pbcshift[i] for i in range(2)]),
						mset=mset_protein,c=clrs[s-1])
		#---inset axis colorbar
		if m == len(analysis_names)-1:
			axins = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation='vertical')
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',
				fontsize=fsaxlabel,rotation=270)
		#---clean up
		if m != 0:
			ax.set_yticklabels([])
			ax.set_ylabel('')
	if len(analysis_names) > 1:
		fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	else: fig.set_size_inches(fig.get_size_inches()[0],fig.get_size_inches()[1])
	gs.tight_layout(fig,h_pad=0.6,w_pad=0.6)
	if len(analysis_names) == 1:
		filespec = specname_guess(sysname,trajsel)
	else:
		filespec = bigname
	plt.savefig(pickles+'fig-topography-'+filespec+'.png',
		dpi=300,bbox_inches='tight')
plt.show()

#---plot the alternative topography correlate curves with more sensible averaging
#---Nb this replaces the original method, which didn't lump frames together
#---Nb the colors and handling the curves for the control case work fine but it is clumsy
if 'topogcorr_plot_alt' in dolist:
	zrange = topogcorr_plot_zrange if len(analysis_names) > 1 else None
	swaprows = True if len(analysis_names) == 1 else False
	fig = plt.figure()
	gs = gridspec.GridSpec((len(analysis_names) if swaprows else 2),
		(2 if swaprows else len(analysis_names)),wspace=0.0,hspace=0.0)
	tcs = []
	axs = []
	avgslices = [slice(0,1),slice(30,60)]
	alpha_slices = [1.,0.2]
	maxyval = 0.
	#---loop over analyses
	for a in analysis_names:
		#---prepare figure with gridspec
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None:
			mset_protein = unpickle(pickles+protein_pkl)
		else:
			mset_protein = mset
		tc = unpickle(pickles+'pkl.topography_correlate.'+specname_guess(sysname,trajsel)+'.pkl')
		tcs.append(tc)
		smoothwid = 4
		ax = fig.add_subplot((gs[0] if len(analysis_names) == 1 
			else gs[(m if swaprows else 0),(0 if swaprows else m)]))
		axs.append(ax)
		ax.set_title(label,fontsize=fsaxtitle)
		#---plot the average structure on the bottom row
		axbot = fig.add_subplot((gs[1] if len(analysis_names) == 1 
			else gs[(m if swaprows else 1),(1 if swaprows else m)]))
		maxval = max(abs(mean(mset.surf,axis=0).min()/10.),abs(mean(mset.surf,axis=0).max()/10.))
		im = plotter2d(axbot,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=cmap,
			inset=False,cmap_washout=1.0,ticklabel_show=[1,1],tickshow=[1,1],
			centertick=False,fs=fsaxlabel,label_style='xy',
			lims=([-maxval,maxval] if zrange == None else [zrange[0],zrange[1]]))
		if m != 0:
			axbot.set_yticklabels([])
			axbot.set_ylabel('')
		#---plot the curves
		hcorr_smooths = []
		for h in range(len(tc.hcorr)):
			#---if not control or one of the names we want to plot
			if (protein_pkl == None and (len(tc.hcorr) == 1 or h > 0)) or \
				(protein_pkl != None and tc.prots_names[h+2] in custom_protein_shifts):
				histdatnorm = nan_to_num([array(tc.topoghist[h][i]/(tc.topoghist[h][i].sum(axis=0))) 
					for i in range(len(tc.topoghist[h]))])
				#---look up the right name for full domain, protein subset, or control
				if protein_pkl == None and h == 0: plotlabel = 'full'
				if protein_pkl == None and h != 0: plotlabel = tc.prots_names[h]
				if protein_pkl != None: plotlabel = tc.prots_names[h+2]
				#---color lookup
				if protein_pkl == None: color = 'k' if h == 0 else clrs[h]
				else:
					if tc.prots_names[h+2] == 'peak': color = clrs[clrdict['red']]
					elif tc.prots_names[h+2] == 'valley': color = clrs[clrdict['blue']]
					else: color = clrs[[i for i in range(len(clrs)) if i not in clrdict.values()][h-2]]
				#---plot curves
				for sl in avgslices:
					dat = mean(sum(histdatnorm,axis=0).T[:,sl],axis=1)
					dat = dat/sum(dat)
					window = float(nhistbins)/nhistbins_smooth
					dat = [sum(dat[i:i+int(window)])/window	for i in range(int(len(dat)-window))]
					maxyval = max(dat) if maxyval < max(dat) else maxyval
					#---suppress the full domain plot for systems with more than one protein
					if (len(tc.hcorr) > 1 and h > 0) or len(tc.hcorr) == 1 or protein_pkl != None:
						ax.plot((array(tc.yedges[:-1]+tc.yedges[1:])/2./10.)[:-window],dat,
							label=plotlabel,lw=3,c=color,alpha=alpha_slices[avgslices.index(sl)])
				#---set protein shifts and plot protein hulls
				if not hasattr(tc,'prots_shifts'): shift = array([0.,0.])
				elif tc.prots_shifts == []: shift = array([0.,0.])
				else: shift = tc.prots_shifts[h]
				for pbcshift in [[0,0],[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]:
					plothull(axbot,mean(mset_protein.protein,axis=0)[tc.prots_refs[h],:2]+\
						shift+array([mean(mset.vecs,axis=0)[i]*pbcshift[i] for i in range(2)]),
						mset=mset_protein,c=color,alpha=1.,ec='w')
		#---inset axis colorbar
		if m == len(analysis_names)-1:
			axins = inset_axes(axbot,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=axbot.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',
				fontsize=fsaxlabel,rotation=270)
		#---clean up
		if m != 0:
			axbot.set_yticklabels([])
			axbot.set_ylabel('')
	for axnum in range(len(axs)):
		ax = axs[axnum]
		ax.set_ylim((0,maxyval*1.1))
		ax.set_xlim((min(tcs[0].yedges)/10.,max(tcs[0].yedges)/10.))
		ax.axvline(x=0,ymin=0.,ymax=1,lw=2,color='k')
		ax.set_xlabel(r'$\left\langle z_r\right\rangle (\mathrm{nm})$',fontsize=fsaxlabel)
		ax.set_yticklabels([])
		ax.grid(True)
		ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=6))			
		x0,x1 = ax.get_xlim()
		y0,y1 = ax.get_ylim()
		ax.set_aspect((x1-x0)/(y1-y0))
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	plt.savefig(pickles+'fig-topography_alt-'+bigname+'.png',dpi=500,bbox_inches='tight')
	if plotviewflag: plt.show()
	plt.clf()
		
#---DEVELOP
#-------------------------------------------------------------------------------------------------------------

#---trying to distinguish control from protein by looking at spread of distributions of bilayer heights
if 0:
	fig = plt.figure()
	for s in range(0,50,5):
		ax = plt.subplot(311)
		ax.plot(sum(sum(tcs[1].topoghist[0][:,s:s+4],axis=0),axis=0),'r',alpha=(1.-s/100.))
	for s in range(0,50,5):
		ax = plt.subplot(312)
		ax.plot(sum(sum(tcs[2].topoghist[1][:,s:s+4],axis=0),axis=0),'b',alpha=(1.-s/100.))
	for s in range(0,50,5):
		ax = plt.subplot(313)
		ax.plot(sum(sum(tcs[2].topoghist[2][:,s:s+4],axis=0),axis=0),'g',alpha=(1.-s/100.))
	plt.show()
if 0:	
	fig = plt.figure()
	ax = plt.subplot(111)
	smoothwind = 2
	ax.plot([std(mean(mean(tcs[1].topoghist[0][:,s:s+smoothwind],axis=0),axis=0)) 
		for s in range(0,70-smoothwind,1)],'ro-')
	ax.plot([std(mean(mean(tcs[2].topoghist[1][:,s:s+smoothwind],axis=0),axis=0)) 
		for s in range(0,70-smoothwind,1)],'bo-')
	ax.plot([std(mean(mean(tcs[2].topoghist[2][:,s:s+smoothwind],axis=0),axis=0)) 
		for s in range(0,70-smoothwind,1)],'go-')
	plt.show()


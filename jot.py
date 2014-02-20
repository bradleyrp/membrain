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

#---plan
analysis_descriptors = {
	'v700-500000-600000-200':
		{'pklfile':
		'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',
		'label':r'$\mathrm{{EXO70}\ensuremath{\times}2{\small (parallel)}}$','nprots':2,
		'custom_topogcorr_specs':[[slice(0,23),'coiled coil'],[slice(1224,1226),
			'PIP2 (1)'],[slice(571,573),'PIP2 (1)']]},
	'v550-300000-400000-200':
		{'pklfile':
		'/structures-broken-transposer-error/pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',
		'label':r'$\mathrm{control}$','nprots':0},
	'v614-120000-220000-200':
		{'pklfile':'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$','nprots':4,
		'topogcorr_pkl':'pkl.topogcorr.membrane-v614-s9-lonestar-120000-220000-200.pkl'}}
analysis_names = ['v700-500000-700000-400','v614-120000-220000-200','v550-300000-400000-200']
analysis_names = ['v614-120000-220000-200']
plot_reord = analysis_names
dolist = ['phase_std_spec2d']
dolist = ['topogcorr_plot']

#---settings
cmap = mpl.cm.RdBu_r

#---settings, topographic correlate
makevid = False #---code needs cleaned up if you use this
nprots = 2
binw = (10,1)

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		msets.append(unpickle(pickles+pklfile))
		msets[-1].calculate_undulations()
		#---hack
		msets[-1].rounder = 20.

#---phase plots
#---Nb this is a canonical spectrum-plotting function, consider standardizing and adding to plotter.py
if 'phase_std_spec2d' in dolist:
	calcs = ['mean','var']
	calcnamesl = ['mean phase','phase variance']
	calcnamesr = [r'$\mathrm{rad}$',r'$(\mathrm{rad}^{2}$)']
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
				ticklabel_show=[(1 if c == len(calcs)-1 else 0),0],tickshow=True,label_style='q')
			#---plot details
			if c == 0: ax.set_title(label)
		#---colorbar for this row
		cbar = gs[c].cbar_axes[0].colorbar(im)
		gs[c].cbar_axes[0].set_ylabel(calcnamesr[c],fontsize=fsaxlabel)
		gs[c][0].set_ylabel(calcnamesl[c],fontsize=fsaxlabel)
	#plt.savefig(pickles+'fig-phases-',dpi=500,bbox_inches='tight')	
	plt.show()

#---DEV
#-------------------------------------------------------------------------------------------------------------

makevid = False #---code needs cleaned up if you use this
nprots = 2
binw = (10,1)

def torus_metric(x1,x2,vecs,ndims=2):
	'''The torus metric i.e. how to measure distance in periodic spaces, quickly.'''
	raw = array([cdist(x1[:,i:i+1],x2[:,i:i+1]) for i in range(ndims)])
	pbcd = array([raw[i]-(raw[i]>vecs[i]/2)*vecs[i] for i in range(ndims)]).T
	return sqrt(sum(pbcd**2,axis=-1)).T

class TopographyCorrelate():
	'''Class object for holding analysis of topography correlation data.'''
	def __init__(self,mset,nprots):
		self.histdat = []
		self.nprots = nprots
		self.prots_nres = shape(mset.protein[0])[0]/nprots
		#---this class is largely just a wrapper for a list of topographic matrices
		self.prots_refs = [slice(None,None)]+[slice(i*self.prots_nres,(i+1)*self.prots_nres) 
			for i in range(nprots)]
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
if 'topogcorr_calc' in dolist and 0:
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		tc = TopographyCorrelate(mset,nprots)
		for i in analysis_descriptors[a]: 
			tc.addnote([i,str((analysis_descriptors[a])[i])])
		for protsref in tc.prots_refs:
			hist_ranges = ((0,ceil(tc.cutoff/100)*100),(-40,40))
			nbins = [(hist_ranges[i][1]-hist_ranges[i][0])/binw[i] for i in range(2)]
			topogmat = []
			topoghist = []
			for fr in range(len(mset.surf)):
				print 'subset: '+str(protsref)+', frame: '+str(fr)
				protpts = mset.protein[fr][protsref,0:2]
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
		special_name = "-".join([i for i in re.match('.+membrane\-.+',
			(analysis_descriptors[a])['pklfile']).string.split('.') 
			if re.match('[a-z][0-9]\-.+',i) 
			or re.match('membrane-v.+',i) 
			or re.match('[0-9]+\-[0-9]+\-[0-9]',i)])
		pickledump(tc,'pkl.topogcorr.'+special_name+'.pkl',directory=pickles)
	
#---method
topogcorr_plot_struct = True

#---DO: add wrapper for subplots here?
axes = []

#---plot the topography correlate curves
if 'topogcorr_plot' in dolist:
	#---prepare figure with gridspec
	fig = plt.figure(figsize=(5*(1+(1 if topogcorr_plot_struct else 0))+3,5))
	gs = gridspec.GridSpec(1,1+(1 if topogcorr_plot_struct else 0))
	#---loop over analyses
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		tc = unpickle(pickles+(analysis_descriptors[a])['topogcorr_pkl'])
		smoothwid = 4
		ax = fig.add_subplot(gs[0])
		for h in range(len(tc.hcorr)):
			#---set labels
			if h == 0:
				label = 'all'
			elif h < tc.nprots+1:
				label = 'domain '+str(h-1)
			else:
				label = None
			#---calculate
			#---Nb since we use the argmax to give the mode of a broad distribution, hard to find error bars
			#---Nb tried taking the second moment relative to the mode but this made no sense
			histdatnorm = nan_to_num([array(tc.topoghist[h][i]/(tc.topoghist[h][i].sum(axis=0))) 
				for i in range(len(tc.topoghist[h]))])
			hcorr = mean(histdatnorm,axis=0)
			hcorr_smooth = [mean([tc.yedges[hcorr[j].argmax()] for j in range(i+smoothwid)]) 
				for i in range(len(hcorr)-smoothwid)]
			ax.plot(tc.xedges[1:-smoothwid],hcorr_smooth,'o-',label=label)
			ax.legend()
			ax.grid(True)
		ax.set_title('Mode(z)',fontsize=fsaxlabel)
		if struct_inset:
			ax = fig.add_subplot(gs[1])
			plotter2d(ax,mset,dat=mean(mset.surf,axis=0),tickshow=True,lognorm=False,cmap=cmap,
				lims=None,inset=False,cmap_washout=0.8,ticklabel_show=False)
			ax.set_title('structure',fontsize=fsaxlabel)
			ax.set_ylabel(r'$y\:(\mathrm{nm})$',fontsize=fsaxlabel)
			ax.set_xlabel(r'$x\:(\mathrm{nm})$',fontsize=fsaxlabel)
	plt.show()

#---DEV

#---plot the 2D histogram
if 0:
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.imshow(hcorr.T,interpolation='nearest',origin='lower',cmap=mpl.cm.binary)
	ax.set_xticks(range(len(xedges))[::4])
	ax.set_yticks(range(len(yedges))[::4])
	ax.set_xticklabels(xedges[::4])
	ax.set_yticklabels(yedges[::4])
	plt.show()

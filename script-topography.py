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

#---plan
analysis_descriptors = {
	'v700-500000-600000-200':
		{'pklfile':'pkl.structures.membrane-v700.u1-lonestar.500000-600000-200.pkl',
		'label':r'$\mathrm{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',
		'nprots':2,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'custom_topogcorr_specs':[[range(0,85)+range(653,653+85),'coiled coil'],[slice(571-1,571-1+2),
			'PIP2 (1)'],[slice(653+571-1,653+571-1+2),'PIP2 (1)']],
		'topogcorr_pkl':'pkl.topogcorr.membrane-v700-u1-lonestar-500000-600000-200.pkl'},
	'v550-400000-500000-160':
		{'pklfile':'pkl.structures.membrane-v550.v1-lonestar.400000-500000-160.pkl',
		'label':r'$\mathrm{control}$',
		'nprots':1,
		'whichframes':slice(0,500),
		'protein_pkl':'pkl.structures.membrane-v612.t4-lonestar.75000-175000-200.pkl',
		'custom_topogcorr_specs':None,
		'custom_protein_shifts':['unshifted','peak','valley'],
		'topogcorr_pkl':'pkl.topogcorr.membrane-v550-v1-lonestar-400000-500000-160.pkl'},
	'v614-120000-220000-200':
		{'pklfile':'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'nprots':4,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'topogcorr_pkl':'pkl.topogcorr.membrane-v614-s9-lonestar-120000-220000-200.pkl'},
	'v612-75000-175000-200':
		{'pklfile':'pkl.structures.membrane-v612.t4-lonestar.75000-175000-200.pkl',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$',
		'nprots':1,
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'whichframes':slice(None,None),
		'topogcorr_pkl':'pkl.topogcorr.membrane-v612-t4-lonestar-75000-175000-200.pkl'}}
do = 'topography'
if do == 'phasecomp':
	analysis_names = ['v700-500000-600000-200','v614-120000-220000-200','v550-400000-500000-160']
	plot_reord = analysis_names
	dolist = ['phase_std_spec2d']
	bigname = '.'.join(analysis_names)
elif do == 'topography':
	analysis_names = ['v614-120000-220000-200','v612-75000-175000-200',
		'v550-400000-500000-160','v700-500000-600000-200']
	plot_reord = analysis_names
	dolist = ['topogcorr_plot']
	bigname = 'v614-v612-v550-v700'

#---settings
cmap = mpl.cm.RdBu_r
clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in range(9)]

#---method, topographic correlate
topogcorr_plot_struct = True

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
		msets.append(unpickle(pickles+pklfile))
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
				centertick=True)
			#---plot details
			if c == 0: ax.set_title(label)
		#---colorbar for this row
		cbar = gs[c].cbar_axes[0].colorbar(im)
		gs[c].cbar_axes[0].set_ylabel(calcnamesr[c],fontsize=fsaxlabel)
		gs[c][0].set_ylabel(calcnamesl[c],fontsize=fsaxlabel)
	plt.savefig(pickles+'fig-phases-'+bigname+'.png',dpi=500,bbox_inches='tight')	
	plt.show()

#---MAIN, PHASE PLOTS
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
		if mset_protein != None:
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
			mset_protein = None
		tc = TopographyCorrelate(mset,nprots,label,mset_protein=mset_protein)
		for i in analysis_descriptors[a]: 
			tc.addnote([i,str((analysis_descriptors[a])[i])])
		#---additional protein sub-domains to test
		if custom_topogcorr_specs != None:
			for item in custom_topogcorr_specs:
				tc.prots_refs.append(item[0])
				tc.prots_names.append(item[1])
		#---if studying a control, get protein points from another pickle and shift them
		elif custom_protein_shifts != None:
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
		for pnum in range(len(tc.prots_refs)):
			protsref = tc.prots_refs[pnum]
			hist_ranges = ((0,ceil(tc.cutoff/100)*100),(-40,40))
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
		pickledump(tc,'pkl.topogcorr.'+specialname((analysis_descriptors[a])['pklfile'])+'.pkl',
			directory=pickles)

#---plot the topography correlate curves
if 'topogcorr_plot_old' in dolist:
	#---loop over analyses
	for a in analysis_names:
		#---prepare figure with gridspec
		fig = plt.figure(figsize=(5*(1+(1 if topogcorr_plot_struct else 0))+3,5))
		gs = gridspec.GridSpec(len(analysis_names),1+(1 if topogcorr_plot_struct else 0))
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None:
			mset_protein = unpickle(pickles+protein_pkl)
		else:
			mset_protein = mset
		tc = unpickle(pickles+(analysis_descriptors[a])['topogcorr_pkl'])
		smoothwid = 4
		ax = fig.add_subplot(gs[0])
		hcorr_smooths = []
		for h in range(len(tc.hcorr)):
			#---set labels
			if h == 0: color = 'k'
			elif h < len(tc.prots_refs): color = clrs[(h-1)%(len(tc.prots_refs)-1)]
			#---calculate
			#---Nb since we use the argmax to give the mode of a broad distribution, hard to find error bars
			#---Nb tried taking the second moment relative to the mode but this made no sense
			histdatnorm = nan_to_num([array(tc.topoghist[h][i]/(tc.topoghist[h][i].sum(axis=0))) 
				for i in range(len(tc.topoghist[h]))])
			hcorr = mean(histdatnorm,axis=0)
			hcorr_smooth = [mean([tc.yedges[hcorr[j].argmax()] for j in range(i+smoothwid)]) 
				for i in range(len(hcorr)-smoothwid)]
			hcorr_smooths.append(hcorr_smooth)
			#---override names if you have a separate protein pickle
			plotlabel = tc.prots_names[(h if protein_pkl == None else h+2)]
			ax.plot(array(tc.xedges[1:-smoothwid])/10.,array(hcorr_smooth)/10.,'o-',label=plotlabel,
				c=color)
			ax.legend()
			ax.grid(True)
		ax.set_xlabel(r'$\min(\left|\mathbf{r}-\mathbf{r}_{\mathrm{prot}}\right|)\:(\mathrm{nm})$',
			fontsize=fsaxlabel)
		ax.set_ylabel(r'$z\:(\mathrm{nm})$',fontsize=fsaxlabel)
		ax.axhline(y=0,xmin=0.,xmax=max(tc.xedges[1:-smoothwid]),lw=2,color='k')
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		if topogcorr_plot_struct:
			ax = fig.add_subplot(gs[1])
			maxval = max(abs(mean(mset.surf,axis=0).min()/10.),abs(mean(mset.surf,axis=0).max()/10.))
			im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=cmap,
				inset=False,cmap_washout=0.8,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-maxval,maxval])
			s0 = 1 if protein_pkl == None else 0
			for s in range((1 if protein_pkl == None else 0),len(tc.prots_refs)):
				color = 'k' if s == 0 else clrs[s-s0]
				for pbcshift in [[0,0],[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]:
					plothull(ax,mean(mset_protein.protein,axis=0)[tc.prots_refs[s],:2]+tc.prots_shifts[s]+\
						array([mean(mset.vecs,axis=0)[i]*pbcshift[i] for i in range(2)]),
						mset=mset_protein,c=clrs[s-1])
			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			fig.colorbar(im,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',fontsize=fsaxlabel)
		plt.setp(cax.get_yticklabels(),fontsize=fsaxlabel)
	plt.savefig(pickles+'fig-topography-'+specialname((analysis_descriptors[a])['pklfile'])+'.png',
		dpi=500,bbox_inches='tight')	
	plt.show()

'''
1. plot multiple systems on one plot
2. plot one plot with multiple structure plots using the right colors
'''

#---plot the topography correlate curves
if 'topogcorr_plot' in dolist:
	#---loop over analyses
	for a in analysis_names:
		#---prepare figure with gridspec
		fig = plt.figure(figsize=(5*(1+(1 if topogcorr_plot_struct else 0))+3,5))
		gs = gridspec.GridSpec(len(analysis_names),1+(1 if topogcorr_plot_struct else 0))
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None:
			mset_protein = unpickle(pickles+protein_pkl)
		else:
			mset_protein = mset
		tc = unpickle(pickles+(analysis_descriptors[a])['topogcorr_pkl'])
		smoothwid = 4
		ax = fig.add_subplot(gs[0])
		hcorr_smooths = []
		for h in range(len(tc.hcorr)):
			#---set labels
			if h == 0: color = 'k'
			elif h < len(tc.prots_refs): color = clrs[(h-1)%(len(tc.prots_refs)-1)]
			#---calculate
			#---Nb since we use the argmax to give the mode of a broad distribution, hard to find error bars
			#---Nb tried taking the second moment relative to the mode but this made no sense
			histdatnorm = nan_to_num([array(tc.topoghist[h][i]/(tc.topoghist[h][i].sum(axis=0))) 
				for i in range(len(tc.topoghist[h]))])
			hcorr = mean(histdatnorm,axis=0)
			hcorr_smooth = [mean([tc.yedges[hcorr[j].argmax()] for j in range(i+smoothwid)]) 
				for i in range(len(hcorr)-smoothwid)]
			hcorr_smooths.append(hcorr_smooth)
			#---override names if you have a separate protein pickle
			plotlabel = tc.prots_names[(h if protein_pkl == None else h+2)]
			ax.plot(array(tc.xedges[1:-smoothwid])/10.,array(hcorr_smooth)/10.,'o-',label=plotlabel,
				c=color)
			ax.legend()
			ax.grid(True)
		ax.set_xlabel(r'$\min(\left|\mathbf{r}-\mathbf{r}_{\mathrm{prot}}\right|)\:(\mathrm{nm})$',
			fontsize=fsaxlabel)
		ax.set_ylabel(r'$z\:(\mathrm{nm})$',fontsize=fsaxlabel)
		ax.axhline(y=0,xmin=0.,xmax=max(tc.xedges[1:-smoothwid]),lw=2,color='k')
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		'''
		if topogcorr_plot_struct:
			ax = fig.add_subplot(gs[1])
			maxval = max(abs(mean(mset.surf,axis=0).min()/10.),abs(mean(mset.surf,axis=0).max()/10.))
			im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=cmap,
				inset=False,cmap_washout=0.8,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-maxval,maxval])
			s0 = 1 if protein_pkl == None else 0
			for s in range((1 if protein_pkl == None else 0),len(tc.prots_refs)):
				color = 'k' if s == 0 else clrs[s-s0]
				for pbcshift in [[0,0],[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]:
					plothull(ax,mean(mset_protein.protein,axis=0)[tc.prots_refs[s],:2]+tc.prots_shifts[s]+\
						array([mean(mset.vecs,axis=0)[i]*pbcshift[i] for i in range(2)]),
						mset=mset_protein,c=clrs[s-1])
			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			fig.colorbar(im,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',fontsize=fsaxlabel)
		plt.setp(cax.get_yticklabels(),fontsize=fsaxlabel)
		'''
	plt.savefig(pickles+'fig-topography-'+specialname((analysis_descriptors[a])['pklfile'])+'.png',
		dpi=500,bbox_inches='tight')	
	plt.show()

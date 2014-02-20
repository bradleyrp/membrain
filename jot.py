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
	'v700-500000-700000-400':
		{'pklfile':
		'/structures-broken-transposer-error/pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',
		'label':r'$\mathrm{{EXO70}\ensuremath{\times}2{\small (parallel)}}$'},
	'v550-300000-400000-200':
		{'pklfile':
		'/structures-broken-transposer-error/pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',
		'label':r'$\mathrm{control}$'},
	'v614-120000-220000-200':
		{'pklfile':'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$'}}
analysis_names = ['v700-500000-700000-400','v614-120000-220000-200','v550-300000-400000-200']
plot_reord = analysis_names
dolist = ['phase_std_spec2d']
dolist = ['test']

#---settings
cmap = mpl.cm.RdBu_r

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
			im = plotter_undulate_spec2d(ax,mset,dat=data,cmap=cmap,lognorm=False,lims=[vmin,vmax],
				ticklabel=[(1 if c == len(calcs)-1 else 0),0],tickshow=True)
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

def torus_metric(x1,x2,vecs,ndims=2):
	'''The torus metric i.e. how to measure distance in periodic spaces, quickly.'''
	raw = array([cdist(x1[:,i:i+1],x2[:,i:i+1]) for i in range(ndims)])
	pbcd = array([raw[i]-(raw[i]>vecs[i]/2)*vecs[i] for i in range(ndims)]).T
	return sqrt(sum(pbcd**2,axis=-1)).T

if 'test' in dolist and 0:
	nbins = 41
	hist_ranges = ((0,400),(-20,20))
	mset = msets[0]
	dat_dist_height = []
	histdat = []
	#---make up 2d histograms
	for fr in range(len(mset.surf)):
		print fr
		protpts = mset.protein[fr][:,0:2]
		surfpts = mset.unzipgrid(mset.surf[fr],vecs=mset.vecs[0])
		#cd = torus_metric(protpts,surfpts,mset.vecs[fr])
		dat_dist_height.append(array([np.min(cd,axis=0),surfpts[:,2]]).T)
		H, xedges, yedges = histogram2d(dat_dist_height[-1][:,0],dat_dist_height[-1][:,1],bins=40,range=hist_ranges)
		row_sums = H.sum(axis=1)
		#histdat.append(nan_to_num(array(H/(H.sum(axis=1))[:, numpy.newaxis])))
		histdat.append(array(H))
	histdat = array(histdat)
	#plotmov(histdat,'test-heightcorr',smooth=1,xedges=[int(i) for i in xedges],yedges=[int(i) for i in yedges])

#plt.imshow(mean(histdat,axis=0).T,interpolation='nearest',origin='lower');plt.show()
if 0:
	sumh = sum(histdat,axis=0)
	sumhn = nan_to_num(sumh/(sumh.sum(axis=0)))
	plt.imshow(sumh.T,interpolation='nearest',origin='lower');plt.show()
donk = [nan_to_num(histdat[i]/histdat[i].sum(axis=1)) for i in range(len(histdat))]
donkm = mean(donk,axis=0)
plt.imshow(donkm.T,interpolation='nearest',origin='lower');plt.show()
#---TOTAL SHITE
#-------------------------------------------------------------------------------------------------------------
'''
if 0:
	dat_dist_height = []
	for fr in range(len(mset.surf)):
		print fr
		protpts = mset.protein[fr][:,0:2]
		surfpts = mset.unzipgrid(mset.surf[fr],vecs=mset.vecs[0])
		cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2])
		dat_dist_height.append(array([np.min(cd,axis=0),surfpts[:,2]]).T)
if 0:
	H, xedges, yedges = histogram2d(mean(dat_dist_height,axis=0)[:,0],
		mean(dat_dist_height,axis=0)[:,1],bins=30,normed=True)
	plt.imshow(array(H).T,interpolation='nearest',origin='lower')
	plt.show()
'''

#!/usr/bin/python -i

#---deprecated by script-postproc-stressmaps.py

#---modified for IET revisions stage 2014.06.13
iet_figs_dir = 'backup-2014.06.12-IET-revision-figs/'
archivedir = './backup-2014.06.01-DEPRECATED/'
structpkldir = 'backup-2014.01.09-enth-review-pickles-exo70-pickles-transpose-error/'

from membrainrunner import *
execfile('locations.py')

import scipy.interpolate
import scipy.integrate
import subprocess

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import ndimage
import matplotlib.gridspec as gridspec

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Nb this script follows script-postproc-stress.py

which_routine = ['exo70pip2_study','exo70pip2_study_v2','enth_study','enth2'][-2]

if which_routine == 'exo70pip2_study':
	#---Note: this is the set for the Exo70+PIP2 simulations
	raw_maps_names = [
			'pkl.stressdecomp.membrane-v701.md.part0003.60000-160000-200.pkl',
			'pkl.stressdecomp.membrane-v700.md.part0002.100000-200000-200.pkl',
			'pkl.stressdecomp.membrane-v550.md.part0006.300000-400000-200.pkl']
	pickle_structure_names = [
			'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
			'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',
			'pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl']
	#---Master ID string
	outname = 'v701.v700.v550.ver3'
	#---Numbers of proteins in each system
	nprots_list = [2,2,0]
	#---Maps labels
	mapslabels = [r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$',
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',
		r'$\textbf{{control}}$']
elif which_routine == 'enth_study':
	#---Note: this is the set for the ENTH simulations
	raw_maps_names = [
		'pkl.stressdecomp.membrane-v614-stress.md.part0002.rerun.pkl',
		'pkl.stressdecomp.membrane-v612-stress.md.part0002.rerun.pkl',
		'pkl.stressdecomp.membrane-v550.md.part0006.300000-400000-200.pkl']
	pickle_structure_names = [
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',
		'pkl.structures.membrane-v612-stress.md.part0003.pkl',
		'pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl']
	outname = 'v614.v612.v550.ver3'	
	#---Numbers of proteins in each system
	nprots_list = [4,1,0]
	#---Maps labels
	mapslabels = [r'$\textbf{{ENTH}\ensuremath{\times}4}$',
		r'$\textbf{{ENTH}\ensuremath{\times}1}$',
		r'$\textbf{{control}}$']
	ncols = len(raw_maps_names)
elif which_routine == 'exo70pip2_study_v2':
	#---Note: this is the set for the Exo70+PIP2 simulations
	raw_maps_names = [
		'pkl.stressdecomp.membrane-v701.md.part0003.60000-160000-200.pkl',
		'pkl.stressdecomp.membrane-v700.md.part0009.500000-700000-400.pkl',
		'pkl.stressdecomp.membrane-v550.md.part0006.300000-400000-200.pkl']
	pickle_structure_names = [
		'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
		'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',
		'pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl']
	#---Master ID string
	outname = 'v701.v700.v550.ver4'
	#---Numbers of proteins in each system
	nprots_list = [2,2,0]
	#---Maps labels
	mapslabels = [r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$',
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',
		r'$\textbf{{control}}$']
elif which_routine == 'enth2':
	#---Note: this is the set for the Exo70+PIP2 simulations
	raw_maps_names = [
		'pkl.stressdecomp.membrane-v701.md.part0003.60000-160000-200.pkl',
		'pkl.stressdecomp.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl',
		'pkl.stressdecomp.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl']
	pickle_structure_names = [
		'pkl.structures.membrane-v700.u1-lonestar-longrun.md.part0009.500000-600000-200.pkl',
		'pkl.structures.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl',
		'pkl.structures.membrane-v550.s0-trajectory-full.md.part0006.300000-400000-200.pkl']
	#---Master ID string
	outname = 'v614.v550'
	#---Numbers of proteins in each system
	nprots_list = [2,4,0]
	#---Maps labels
	mapslabels = [r'$\textbf{{v701trash}}$',r'$\textbf{{ENTH}\ensuremath{\times}4}$',
		r'$\textbf{{control}}$']
	ncols = len(raw_maps_names)

#---plots
plot_maps = 1
plot_hist = 0 ######## needs centered at 0!?!
plot_hist_subdivide = 0
plot_hist_subdivide_mean = 0

#---methods
smoothmaps = False
smoothwindow = 2.

#---settings
nbins = 31

#---file names
plot_maps_file = 'fig-'+outname+'-stress-framewise-maps.png'
plot_hist_file = 'fig-'+outname+'-stress-framewise-histograms.png'

#---sign change
signchange = 1

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'raw_maps' not in globals():
	raw_maps = []
	for name in raw_maps_names:
		print 'Loading maps from '+name
		raw_maps.append(pickle.load(open(pickles+archivedir+name,'r')))
	msets = []
	for name in pickle_structure_names:
		print 'Loading structures from '+name
		msets.append(pickle.load(open(pickles+structpkldir+name,'r')))

#---plot spontaneous curvature (C0) histograms, averaged across all voxels and all frames (together)
if plot_hist:
	which_brewer_colors = [0,1,2,3,4,5,6,7]
	clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
	pdist0 = signchange*array([i for j in mean([i[0] for i in raw_maps[0]],axis=0) for i in j])
	pdist1 = signchange*array([i for j in mean([i[0] for i in raw_maps[1]],axis=0) for i in j])
	pdist2 = signchange*array([i for j in mean([i[0] for i in raw_maps[2]],axis=0) for i in j])
	maxval = max([max(pdist0),max(pdist1),max(pdist2)])
	minval = min([min(pdist0),min(pdist1),min(pdist2)])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.rc('font', family='sans-serif')
	ax.grid(True)
	nbins = 50
	hist0,binedge0 = numpy.histogram(pdist0,bins=nbins,normed=True,range=(minval,maxval))
	hist1,binedge1 = numpy.histogram(pdist1,bins=nbins,normed=True,range=(minval,maxval))
	hist2,binedge2 = numpy.histogram(pdist2,bins=nbins,normed=True,range=(minval,maxval))
	mid0 = (binedge0[1:]+binedge0[:-1])/2
	mid1 = (binedge1[1:]+binedge1[:-1])/2
	mid2 = (binedge2[1:]+binedge2[:-1])/2	
	plt.plot(mid0,hist0,'o-',c=clrs[1],alpha=1.,lw=2,label=mapslabels[0])
	plt.plot(mid1,hist1,'o-',c=clrs[3],alpha=1.,lw=2,label=mapslabels[1])
	plt.plot(mid2,hist2,'o-',c=clrs[5],alpha=1.,lw=2,label=mapslabels[2])
	plt.xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',labelpad = 10,fontsize=20)
	plt.ylabel('frequency', labelpad = 10,fontsize=20)
	ax.set_xlim((-0.05,0.05))
	fig.tight_layout()
	plt.legend()
	plt.tight_layout() 
	plt.savefig(pickles+iet_figs_dir+plot_hist_file,dpi=500,bbox_inches='tight')
	plt.show()
	plt.cla()
	plt.close()
	
#---plot spontaneous curvature (C0) histograms on a selection of the box
if plot_hist_subdivide:
	#---zoom full vs center box, with possible fold-over, all in one plot
	if 0:
		which_brewer_colors = [1,3,5,0,2,4]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))],
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))]]
		zoomnames = (', full',', center')
		fig = plt.figure(figsize=(6,4))
		gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
		lims = (-0.1,0.1)
		plotlims = (-0.06,0.06)
		plotlims = (-.15,.15)
		clrsi = 0
		ax = fig.add_subplot(gs[0])
		for z in range(len(zooms)):
			zoom = zooms[z]
			#for d in range(len(raw_maps)):
			for d in range(0,1):
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-.',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				ax.plot(mid,hist,'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])				
				clrsi += 1
				ax.legend()
			ax.set_xlim(plotlims)
		plt.show()
	#---stacked method without averaging
	if 0:
		plot_hist_subdivide_mean = False
		which_brewer_colors = [1,3,5,0,2,4]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))],
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))]]
		zoomnames = (', full',', center')
		fig = plt.figure(figsize=(6,4))
		gs = gridspec.GridSpec(3,1,wspace=0.0,hspace=0.0)
		lims = (-1,1)
		plotlims = (-0.06,0.06)
		plotlims = (-1,1)
		clrsi = 0
		maxval = 0
		axes = []
		for d in range(len(raw_maps)):
			ax = fig.add_subplot(gs[d])
			axes.append(ax)
			ax.axvline(x=0,ls='-',lw=1,c='k')
			ax.set_ylim((0,0.04))
			for z in range(len(zooms)):
				zoom = zooms[z]
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				print mean(pdist)
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-.',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				ax.plot(mid,hist,'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])				
				clrsi += 1
				ax.legend()
				if max(hist) > maxval: maxval = max(hist)
			ax.set_xlim(plotlims)
		for ax in axes:
			ax.set_ylim((0,maxval*1.1))
		plt.show()
	#---stacked method without averaging next to unstacked method with averaging
	if 0:
		plot_hist_subdivide_mean = False
		which_brewer_colors = [0,1,2,3,4,5,6,7]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))],
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))]]
		zoomnames = (', full',', center')
		fig = plt.figure(figsize=(16,10))
		gs = gridspec.GridSpec(3,2,wspace=0.0,hspace=0.0)
		lims = (-1,1)
		plotlims = (-0.5,0.5)
		clrsi = 0
		maxval = 0
		axes1 = []
		for d in range(len(raw_maps)):
			ax = fig.add_subplot(gs[d,0])
			axes1.append(ax)
			ax.axvline(x=0,ls='-',lw=1,c='k')
			ax.set_ylim((0,0.04))
			for z in range(len(zooms)):
				zoom = zooms[z]
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				print mean(pdist)
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-.',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				ax.plot(mid,hist,'o-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])				
				clrsi += 1
				ax.legend(prop={'size':14})
				if max(hist) > maxval: maxval = max(hist)
			ax.set_xlim(plotlims)
		axes1[0].set_title('all $\mathsf{C_{0}\,(nm^{-1})}$')
		for ax in axes1:
			ax.set_ylim((0,maxval*1.1))
			ax.set_ylabel('frequency',fontsize=18)
			ax.grid(True)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=6,prune='both'))
		axes1[-1].set_xlabel('$\mathsf{C_{0}\,(nm^{-1})}$',fontsize=18)
		plot_hist_subdivide_mean = True
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))],
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))]]
		zoomnames = (', full',', center')
		lims = (-0.1,0.1)
		plotlims = (-0.06,0.06)
		plotlims = (-0.07,0.07)
		clrsi = 0
		maxval = 0
		axes2 = []
		for d in range(len(raw_maps)):
			ax = fig.add_subplot(gs[d,1])
			axes2.append(ax)
			ax.axvline(x=0,ls='-',lw=1,c='k')
			ax.set_ylim((0,0.04))
			for z in range(len(zooms)):
				zoom = zooms[z]
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				print mean(pdist)
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-.',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				ax.plot(mid,hist,'o-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])				
				clrsi += 1
				ax.legend(prop={'size':14})
				if max(hist) > maxval: maxval = max(hist)
			ax.set_xlim(plotlims)
		axes2[0].set_title(r'time-averaged $\mathsf{C_{0}\,(nm^{-1})}$')
		for ax in axes2:
			ax.set_ylim((0,maxval*1.1))
			ax.yaxis.tick_right()
			ax.set_ylabel('frequency',fontsize=18)
			ax.yaxis.set_label_position("right")
			ax.grid(True)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=6,prune='both'))
		axes2[-1].set_xlabel('$\mathsf{C_{0}\,(nm^{-1})}$',fontsize=18)
		plt.savefig(pickles+iet_figs_dir+'fig-stress-master-summary-centerbox-ENTH.png',dpi=500,bbox_inches='tight')
		plt.show()
	#---stacked method WITHOUT averaging WITHOUT THE other one
	if 1:
		which_brewer_colors = [0,1,2,3,4,5,6,7]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		plot_hist_subdivide_mean = True
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))],
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))]]
		zoomnames = (', full',', center')
		lims = (-0.1,0.1)
		plotlims = (-0.06,0.06)
		clrsi = 0
		maxval = 0
		axes2 = []
		fig = plt.figure(figsize=(6,8))
		gs = gridspec.GridSpec(3,1,wspace=0.0,hspace=0.0)
		for d in range(len(raw_maps)):
			ax = fig.add_subplot(gs[d,0])
			axes2.append(ax)
			ax.axvline(x=0,ls='-',lw=1,c='k')
			ax.set_ylim((0,0.04))
			for z in range(len(zooms)):
				zoom = zooms[z]
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				print mean(pdist)
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-.',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])
				ax.plot(mid,hist,'o-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+zoomnames[z])				
				clrsi += 1
				ax.legend(prop={'size':14})
				if max(hist) > maxval: maxval = max(hist)
			ax.set_xlim(plotlims)
		for ax in axes2:
			ax.set_ylim((0,maxval*1.1))
			#ax.yaxis.tick_right()
			ax.set_ylabel('frequency',fontsize=18)
			#ax.yaxis.set_label_position("right")
			ax.grid(True)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=7,prune='both'))
		for ax in axes2[:-1]:
			ax.set_xticklabels([])
		axes2[-1].set_xlabel('$\mathsf{C_{0}\,(nm^{-1})}$',fontsize=18)
		plt.savefig(pickles+iet_figs_dir +'fig-stress-master-summary-centerbox-ENTH.png',dpi=500,bbox_inches='tight')
		plt.show()
	#---new method breaks it down into (1) positive vs negative (2) positive center vs full (3) center v full
	if 0:
		which_brewer_colors = [0,2,4,1,3,5]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))],
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))]]
		zoomnames = (', center',', full')
		fig = plt.figure(figsize=(6,4))
		gs = gridspec.GridSpec(3,2,wspace=0.0,hspace=0.0)
		lims = (-0.1,0.1)
		#---full zoom, positive vs negative
		ax = fig.add_subplot(gs[0,0:2])
		z = 0
		zoom = zooms[z]
		clrsi = 0
		for d in range(len(raw_maps)):
			dat = raw_maps[d]
			if plot_hist_subdivide_mean:
				pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
			else:
				pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
			hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
				range=lims)
			mid = (binedge[1:]+binedge[:-1])/2
			posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
			negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
			ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			clrsi += 1
			ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			clrsi += 1
			ax.legend()
		#---now positive full vs center
		clrsi = 0
		ax = fig.add_subplot(gs[1,0:2])
		for z in range(len(zooms)):
			zoom = zooms[z]
			for d in range(len(raw_maps)):
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				#negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				clrsi += 1
				ax.legend()
			ax.set_xlim(plotlims)
		#---now negative full vs center
		clrsi = 0
		ax = fig.add_subplot(gs[2,0:2])
		for z in range(len(zooms)):
			zoom = zooms[z]
			for d in range(len(raw_maps)):
				dat = raw_maps[d]
				if plot_hist_subdivide_mean:
					pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				else:
					pdist = signchange*flatten(array([array(i[0])[zoom[0],zoom[1]] for i in dat]))
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				#posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c=clrs[clrsi%len(clrs)],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				clrsi += 1
				ax.legend()
			ax.set_xlim(plotlims)
		plt.show()
	#---new method does the same thing, but also breaks it out by system
	if 0:
		which_brewer_colors = [0,2,4,1,3,5]
		clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		gridsize = shape(raw_maps[0][0][0])[0]
		zooms = [
			[slice(int(round(0.*gridsize)),int(round(1.*gridsize))),
				slice(int(round(0.*gridsize)),int(round(1.*gridsize)))],
			[slice(int(round(0.25*gridsize)),int(round(0.75*gridsize))),
				slice(int(round(0.25*gridsize)),int(round(0.75*gridsize)))]]
		zoomnames = (', center',', full')
		fig = plt.figure(figsize=(6,4))
		gs = gridspec.GridSpec(3,3,wspace=0.0,hspace=0.0)
		lims = (-0.1,0.1)
		plotlims = (0,0.04)
		z = 0
		zoom = zooms[0]
		#---top row is full system positive vs negative
		for col in range(ncols):
			dat = raw_maps[col]
			ax = fig.add_subplot(gs[0,col])
			pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
			hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
				range=lims)
			mid = (binedge[1:]+binedge[:-1])/2
			posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
			negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
			ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c='r',alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c='b',alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			ax.set_xlim(plotlims)
			ax.set_ylim((0,0.4))		
		#---middle row is positive curvatures
		clrsi = 0
		for col in range(ncols):
			dat = raw_maps[col]
			ax = fig.add_subplot(gs[1,col])
			for z in range(len(zooms)):
				zoom = zooms[z]
				pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c=clrs[z+4],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				#ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c='b',alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			ax.set_xlim(plotlims)
			ax.set_ylim((0,0.4))
		#---bottom row
		clrsi = 0
		for col in range(ncols):
			dat = raw_maps[col]
			ax = fig.add_subplot(gs[2,col])
			for z in range(len(zooms)):
				zoom = zooms[z]
				pdist = signchange*flatten(mean(array([array(i[0])[zoom[0],zoom[1]] for i in dat]),axis=0))
				hist,binedge = numpy.histogram(pdist,bins=nbins,weights=[1./len(pdist) for i in pdist],
					range=lims)
				mid = (binedge[1:]+binedge[:-1])/2
				posmid = [i for i in range(len(hist)) if round(mid[i],8) >= 0.]
				negmid = [i for i in range(len(hist)) if round(mid[i],8) <= 0.]
				#ax.plot([mid[i] for i in posmid],[hist[i] for i in posmid],'-',c='r',alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
				ax.plot([-mid[i] for i in negmid],[hist[i] for i in negmid],'-',c=clrs[z+4],alpha=1.,lw=2,label=mapslabels[d]+' '+zoomnames[z])
			ax.set_xlim(plotlims)
			ax.set_ylim((0,0.4))		
		plt.show()




#---plot spontaneous curvature (C0) maps, averaged over all frames
if plot_maps:
	prot_centers = []
	for m in range(len(msets)):
		mset = msets[m]
		nprots = nprots_list[m]
		if nprots > 0:
			protlen = int(shape(mset.protein[0])[0]/nprots)
			protlenshow = protlen
			prot_centers.append(mean([mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow]
				for i in range(nprots)],axis=1))
		else:
			prot_centers.append(None)
	vecs = mean(msets[0].vecs,axis=0)
	span,nnum,numgridpts,distance_factor = [4,32,64,1]
	result_stack = [[mean([i[0] for i in raw_maps[j]],axis=0),prot_centers[j]] for j in range(len(raw_maps))]
	fig = plt.figure()	
	gs = mpl.gridspec.GridSpec(4,1,width_ratios=[1,1,2,1],height_ratios=[1])
	plt.rc('font', family='sans-serif')
	extremum = max([max([max(i) for i in result_stack[j][0]]) for j in range(ncols)])
	extremum = 0.05
	ax0 = plt.subplot2grid((1,4),(0,0))
	ax0.set_title(mapslabels[0])
	ax0.set_xticklabels([])
	ax0.set_yticklabels([])
	ax0.set_adjustable('box-forced')

	dat = result_stack[0][0]
	if smoothmaps:
		dat2 = []
		for i in dat:
			dat2.append(ndimage.gaussian_filter(i, sigma=smoothwindow*256/(4.*50)))
		dat = dat2
	
	nprots = nprots_list[0]
	if nprots > 0:
		protlen = int(shape(msets[0].protein[0])[0]/nprots)
		for protpts in [mean(msets[0].protein,axis=0)[i*protlen:i*protlen+protlenshow] 
			for i in range(nprots)]:
			protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in protpts])
			hull = scipy.spatial.ConvexHull(protpts)
			ax0.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'k-',lw=0.6)
			shifthullx = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][0] 
				for i in range(len(hull.vertices))]
			shifthully = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][1] 
				for i in range(len(hull.vertices))]
			ax0.plot(shifthullx,shifthully,'k-',lw=1.0)
	ax0.imshow(signchange*array(dat).T,interpolation='nearest',
		origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax1 = plt.subplot2grid((1,4),(0,1))
	ax1.set_title(mapslabels[1])
	ax1.set_xticklabels([])
	ax1.set_yticklabels([])
	ax1.set_adjustable('box-forced')

	dat = result_stack[1][0]
	if smoothmaps:
		dat2 = []
		for i in dat:
			dat2.append(ndimage.gaussian_filter(i, sigma=smoothwindow*256/(4.*50)))
		dat = dat2

	nprots = nprots_list[1]
	if nprots > 0:
		protlen = int(shape(msets[1].protein[0])[0]/nprots)
		for protpts in [mean(msets[1].protein,axis=0)[i*protlen:i*protlen+protlenshow] 
			for i in range(nprots)]:
			protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in protpts])
			hull = scipy.spatial.ConvexHull(protpts)
			ax1.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'k-',lw=0.6)
			shifthullx = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][0] 
				for i in range(len(hull.vertices))]
			shifthully = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][1] 
				for i in range(len(hull.vertices))]
			ax1.plot(shifthullx,shifthully,'k-',lw=1.0)
	ax1.imshow(signchange*array(dat).T,interpolation='nearest',origin='LowerLeft',
		vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	
	ax2 = plt.subplot2grid((1,4), (0,2),colspan=1)
	ax2.set_title(mapslabels[2])
	ax2.set_xticklabels([])
	ax2.set_yticklabels([])
	ax2.set_adjustable('box-forced')
	dat = result_stack[2][0]
	
	if smoothmaps:
		dat2 = []
		for i in dat:
			dat2.append(ndimage.gaussian_filter(i, sigma=smoothwindow*256/(4.*50)))
		dat = dat2

	img = ax2.imshow(-1*array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	cax = inset_axes(ax2,
		width="5%",
		height="100%",
		bbox_transform=ax2.transAxes,
		bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
		loc= 1)
	fig.colorbar(img,cax=cax)
	cax.tick_params(labelsize=10) 
	cax.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=10)
	plt.tight_layout()
	plt.savefig(pickles+iet_figs_dir+plot_maps_file,dpi=500,bbox_inches='tight')
	plt.show()
	plt.cla()
	plt.close()


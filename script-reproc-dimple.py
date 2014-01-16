#!/usr/bin/python -i

from membrainrunner import *

import os
from scipy.optimize import leastsq
import matplotlib as mpl
from pylab import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

skip = 1
framecount = None
location = ''
execfile('locations.py')

import matplotlib.gridspec as gridspec
which_brewer_colors = [0,1,2,3,4,5,6,7]
clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
clrs2 = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams.update({'font.size': 14})

analysis_descriptors = [
	('pkl.dimple.v614-stress.md.part0002.rerun.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4}$',1,
		'pkl.tilefilter-areas.v614-stress.md.part0002.rerun.pkl'),
	('pkl.dimple.v612-stress.md.part0003.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{ENTH}\ensuremath{\times}1}$',1,
		'pkl.tilefilter-areas.v612-stress.md.part0003.pkl'),
	('pkl.dimple.v612-stress.md.part0003.prot-v614.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{ENTH}\ensuremath{\times}1(4x area)}$',1,
		'pkl.tilefilter-areas.v612-stress.md.part0003.prot-v614.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),

	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.invert.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control, invert}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),


	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v700.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v700.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),


	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v700.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v700.pkl'),
	('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',1,
		'pkl.tilefilter-areas.v700.md.part0002.100000-200000-200.pkl'),
	('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$',1,
		'pkl.tilefilter-areas.v701.md.part0003.60000-160000-200.pkl')]

plotspecs = 'enth'
plotspecs = 'both'
if plotspecs == 'both':
	analysis_plan = slice(None,None)
	appor = (0,1,2,3,4,5,5,5,5,6,7,8)
	figoutname = 'fig-dimple-master-summary-ENTH-EXO70.png'
	figsize = (14,16)
elif plotspecs == 'enth':
	analysis_plan = slice(0,3)
	appor = (0,1,2)
	figoutname = 'fig-dimple-master-summary-ENTH.png'
	figsize = (14,8)

do_stacked_plot = True
do_stacked_plot_with_sigma = True
do_stacked_plot_ver1 = False
do_opposite_signs = False
do_single_plot = False

results_stack = []
for pnum in range(len(analysis_descriptors[analysis_plan])):
	results_stack.append(pickle.load(open(pickles+analysis_descriptors[analysis_plan][pnum][0],'r')))
results_areas_stack = []
for pnum in range(len(analysis_descriptors[analysis_plan])):
	results_areas_stack.append(pickle.load(open(pickles+analysis_descriptors[analysis_plan][pnum][4],'r')))

	
nbins = 20
nbins_sigma = 20
minval,maxval = -0.10,0.10
minval_sigma,maxval_sigma = 0,30

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Advanced plotting method
if do_stacked_plot:
	fig = plt.figure(figsize=figsize)
	gs = gridspec.GridSpec(max(appor)+1,5,wspace=0.0,hspace=0.0)
	#---plot maximum mean curvatures	
	maxpeak = 0
	axes_maxcurv = []
	for p in range(len(analysis_descriptors[analysis_plan])):
		if appor[p] > 0 and appor[p] == appor[p-1]:
			thisaxis = axes_maxcurv[-1]
		else:
			thisaxis = fig.add_subplot(gs[appor[p],0:2])
		ccodes = analysis_descriptors[analysis_plan][p][1]
		name = analysis_descriptors[analysis_plan][p][2]
		fillcode = analysis_descriptors[analysis_plan][p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		for o in range(len(order)):
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=False,
				weights=[1./len(validhs) for i in validhs],range=(minval,maxval))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			if o == 1:
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1.,lw=2,label=name)
			elif o == 0:
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1.,lw=2)
			else:
				thisaxis.plot(mid0,hist0,'--',c='k',alpha=1.,lw=2)
			if fillcode and o != 2:
				thisaxis.fill_between(mid0,hist0,[0 for i in mid0],facecolor=ccodes[o],
					alpha=0.2,interpolate=True)
			if max(hist0) > maxpeak: maxpeak = max(hist0)
		if not (appor[p] > 0 and appor[p] == appor[p-1]):
			axes_maxcurv.append(thisaxis)
	for a in range(len(axes_maxcurv)):
		ax = axes_maxcurv[a]
		if a == 0:
			ax.set_title(r'$\textbf{mean curvatures}$')
		ax.set_ylim(0,1.2*maxpeak)
		ax.axvline(x=0,ls='-',lw=1,c='k')
		ax.legend(loc=2,prop={'size':10})
		ax.set_ylabel('frequency')
		ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
		ax.grid(True)
		if a == len(axes_maxcurv)-1:
			ax.spines['bottom'].set_position(('outward', 10))
			ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
			second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
			ax.spines['second_bottom'] = second_bottom
			ax.set_xticks(arange(-0.10,0.12,0.04))
		else:
			ax.set_xticks(arange(-0.10,0.12,0.04))
			ax.set_xticklabels([])
	#---plot extents of curvature
	axes_sigmas = []
	maxpeak = 0
	for p in range(len(analysis_descriptors[analysis_plan])):
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = [1] if expected_direction == 1 else [0]
		o = 0
		thisaxis = plt.subplot(gs[appor[p],2])
		params = results_stack[p][order[o]].get(['type','params'])
		maxhs = results_stack[p][order[o]].get(['type','maxhs'])
		ccodes = analysis_descriptors[analysis_plan][p][1]
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 
			and abs(10*maxhs[i]) < 0.1)]
		sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0]
		sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0]
		hist0,binedge0 = numpy.histogram(sigma_x,bins=nbins_sigma,normed=True,density=True,
			range=(minval_sigma,maxval_sigma))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		thisaxis.plot(mid0,hist0,c=ccodes[1],alpha=1.,lw=2)
		hist1,binedge1 = numpy.histogram(sigma_y,bins=nbins_sigma,normed=True,density=True,
			range=(minval_sigma,maxval_sigma))
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		thisaxis.plot(mid1,hist1,c=ccodes[0],alpha=1.,lw=2)
		axes_sigmas.append(thisaxis)
		if max(max(hist0),max(hist1)) > maxpeak: maxpeak = max(max(hist0),max(hist1))
	for a in range(len(axes_sigmas)):
		ax = axes_sigmas[a]
		if a == 0:
			ax.set_title(r'$\textbf{extents}$')
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
		ax.set_xticks(arange(5,30,5))
		if a == len(axes_sigmas)-1:
			ax.set_xlabel('$\mathsf{\sigma_a,\sigma_b\,(nm^{2})}$',fontsize=14)
		else:
			ax.set_xticklabels([])		
	#---plot areas
	axes_areas = []
	maxpeak = 0
	for p in range(len(analysis_descriptors[analysis_plan])):
		if appor[p] > 0 and appor[p] == appor[p-1]:
			thisaxis = axes_areas[-1]
			repeat = True
		else:
			thisaxis = fig.add_subplot(gs[appor[p],3])
			repeat = False
		area_per_tile = results_areas_stack[p].notes[[i[0] 
			for i in results_areas_stack[p].notes].index('area_per_tile')][1]
		area_counts = results_areas_stack[p].data
		posarea = array([area_per_tile*area_counts[i][2] for i in range(len(area_counts))])
		negarea = array([area_per_tile*area_counts[i][3] for i in range(len(area_counts))])
		thisaxis.plot(posarea,'r-',label=(None if repeat else '$z>0$' ),lw=1,alpha=0.8)
		thisaxis.plot(negarea,'b-',label=(None if repeat else '$z<0$'),lw=1,alpha=0.8)
		t = range(len(area_counts))
		thisaxis.fill_between(t, negarea,posarea, facecolor='b',alpha=0.25,where=negarea>posarea)
		thisaxis.fill_between(t, posarea,negarea, facecolor='r',alpha=0.25,where=negarea<posarea)
		axes_areas.append(thisaxis)
		if max(max(posarea),max(negarea)) > maxpeak: maxpeak = max(max(posarea),max(negarea))
	for a in range(len(axes_areas)):
		ax = axes_areas[a]
		if a == 0:
			ax.set_title(r'$\textbf{areas (+/-)}$')
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=6))
		if a == len(axes_areas)-1:
			ax.set_xlabel('frame',fontsize=14)
			outlegend = ax.legend(bbox_to_anchor=(0.5,-0.6),
				loc=8,borderaxespad=-1.,prop={'size':10},bbox_transform=ax.transAxes)
		else:
			ax.set_xticklabels([])		
	#---plot area histograms
	axes_area_hists = []
	minval_areas = 0
	maxval_areas = maxpeak
	maxfreq = 0
	for p in range(len(analysis_descriptors[analysis_plan])):
		if appor[p] > 0 and appor[p] == appor[p-1]:
			thisaxis = axes_area_hists[-1]
		else:
			thisaxis = fig.add_subplot(gs[appor[p],4])
		area_per_tile = results_areas_stack[p].notes[[i[0] 
			for i in results_areas_stack[p].notes].index('area_per_tile')][1]
		area_counts = results_areas_stack[p].data
		posarea = array([area_per_tile*area_counts[i][2] for i in range(len(area_counts))])
		negarea = array([area_per_tile*area_counts[i][3] for i in range(len(area_counts))])
		hist0,binedge0 = numpy.histogram(posarea,bins=nbins_sigma,normed=False,
			weights=[1./len(negarea) for i in negarea],
			range=(minval_areas,maxval_areas))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		thisaxis.plot(hist0,mid0,c='r',alpha=1.,lw=1)
		hist1,binedge1 = numpy.histogram(negarea,bins=nbins_sigma,normed=False,
			weights=[1./len(negarea) for i in negarea],
			range=(minval_areas,maxval_areas))
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		thisaxis.plot(hist1,mid1,c='b',alpha=1.,lw=1)
		thisaxis.fill_betweenx(mid0,[0 for i in mid0],hist0,facecolor='r',
			alpha=0.2)
		thisaxis.fill_betweenx(mid1,[0 for i in mid1],hist1,facecolor='b',
			alpha=0.2)
		axes_area_hists.append(thisaxis)
		if max(max(hist0),max(hist1)) > maxfreq: maxfreq = max(max(hist0),max(hist1))
	for a in range(len(axes_area_hists)):
		ax = axes_area_hists[a]
		if a == 0:
			ax.set_title(r'$\textbf{areas (+/-)}$')
		ax.grid(True)
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		ax.set_ylabel(r'$\textbf{area (nm\ensuremath{{}^{2}})}$',fontsize=12)
		maxfreq = 0.5
		ax.set_xlim(0,1.0*maxfreq)
		ax.get_yaxis().set_major_locator(MaxNLocator(prune='both',nbins=6))
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=6))
		if a == len(axes_area_hists)-1:
			ax.set_xlabel('frequency',fontsize=14)
		else:
			ax.set_xticklabels([])
	plt.subplots_adjust(hspace = 0)
	plt.subplots_adjust(wspace = 0)
	plt.savefig(pickles+figoutname,dpi=300,bbox_extra_artists=(outlegend,), bbox_inches='tight')
	#---print report
	for p in range(len(analysis_descriptors[analysis_plan])):
		ccodes = analysis_descriptors[analysis_plan][p][1]
		name = analysis_descriptors[analysis_plan][p][2]
		print name
		fillcode = analysis_descriptors[analysis_plan][p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		for o in range(len(order)):
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0]
			print mean(validhs)
			print mean(sigma_x)
			print mean(sigma_y)
	plt.show()	


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

#---deprecated code for the original plotting method
'''
analysis_descriptors = [
	('pkl.dimple.v614-stress.md.part0002.rerun.pkl',),
	('pkl.dimple.v612-stress.md.part0003.pkl',),
	('pkl.dimple.v550.md.part0006.300000-400000-200.pkl',),
	('pkl.dimple.v550.md.part0006.300000-400000-200.testshift-v614.pkl',),
	('pkl.dimple.v550.md.part0006.300000-400000-200.testshift-v700.pkl',),
	('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',),
	('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',)]
analysis_plan = range(len(analysis_descriptors))
names = (r'$\textbf{{ENTH}\ensuremath{\times}4}$',r'$\textbf{{ENTH}\ensuremath{\times}1}$',
	r'$\textbf{control(1)}$',r'$\textbf{control(2)}$',r'$\textbf{control(3)}$',
	r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$',r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$')
appor = (0,1,2,2,2,3,4)
ccodes = [(clrs[0],clrs[1]),(clrs[2],clrs[3]),
	(clrs[4],clrs[5]),(clrs[6],clrs[7]),(clrs[0],clrs[1]),(clrs[2],clrs[3]),
	(clrs[0],clrs[1]),(clrs[2],clrs[3])]
fillcodes = (1,1,0,0,0,1,1)
'''

analysis_descriptors = [
	('pkl.dimple.v614-stress.md.part0002.rerun.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4}$',1,
		'pkl.tilefilter-areas.v614-stress.md.part0002.rerun.pkl'),
	('pkl.dimple.v612-stress.md.part0003.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{ENTH}\ensuremath{\times}1}$',1,
		'pkl.tilefilter-areas.v612-stress.md.part0003.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control(1)}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v614.pkl'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v700.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control(2)}$',0,
		'pkl.tilefilter-areas.v550.md.part0006.300000-400000-200.prot-v700.pkl'),
	('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',1,
		'pkl.tilefilter-areas.v700.md.part0002.100000-200000-200.pkl'),
	('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$',1,
		'pkl.tilefilter-areas.v701.md.part0003.60000-160000-200.pkl')]
analysis_plan = slice(None,None)
appor = (0,1,2,2,3,4)

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

#---Original plotting method
if do_stacked_plot_ver1:
	fig, axes = plt.subplots(nrows=max(appor)+1,ncols=(2 if do_stacked_plot_with_sigma else 1),
		figsize=((15 if do_stacked_plot_with_sigma else 10),2*(max(appor)+1)))
	maxpeak = 0
	for p in range(len(analysis_descriptors[analysis_plan])):
		thisaxis = axes[appor[p]][0] if do_stacked_plot_with_sigma else axes[appor[p]]
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
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
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
	for p in range(max(appor)+2):
		if p == 1:
			thisaxis.set_title(r'$\textbf{mean curvatures}$')
		thisaxis = axes[appor[p]][0] if do_stacked_plot_with_sigma else axes[appor[p]]
		thisaxis.set_ylim(0,1.1*maxpeak)
		thisaxis.axvline(x=0,ls='-',lw=1,c='k')
		thisaxis.legend(loc=2,fontsize=12)
		thisaxis.set_ylabel('frames')
		thisaxis.get_yaxis().set_major_locator(MaxNLocator(prune='both'))
		thisaxis.grid(True)
		thisaxis.locator_params(nbins=4)
		if p == max(appor)+1:
			thisaxis.set_xlabel('$\mathsf{H_{max}(nm^{-1})}$',fontsize=14)
		else:
			thisaxis.set_xticklabels([])
	if do_stacked_plot_with_sigma:	
		for p in range(len(analysis_descriptors[analysis_plan])):
			expected_direction = results_stack[p][0].notes[([i[0] 
				for i in results_stack[p][0].notes].index('expected_direction'))][1]
			order = [1] if expected_direction == 1 else [0]
			o = 0
			thisaxis = axes[appor[p]][1]
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			ccodes = analysis_descriptors[analysis_plan][p][1]
			validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 
				and abs(10*maxhs[i]) < 0.1)]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0]
			hist0,binedge0 = numpy.histogram(sigma_x,bins=nbins_sigma,normed=True,
				range=(minval_sigma,maxval_sigma))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			thisaxis.plot(mid0,hist0,c=ccodes[1],alpha=1.,lw=2)
			hist1,binedge1 = numpy.histogram(sigma_y,bins=nbins_sigma,normed=True,
				range=(minval_sigma,maxval_sigma))
			mid1 = (binedge1[1:]+binedge1[:-1])/2
			thisaxis.plot(mid1,hist1,c=ccodes[0],alpha=1.,lw=2)
		for p in range(max(appor)+2):
			if p == 1:
				thisaxis.set_title(r'$\textbf{extents of curvature}$')
			thisaxis = axes[appor[p]][1] if do_stacked_plot_with_sigma else axes[appor[p]]
			thisaxis.grid(True)
			thisaxis.yaxis.tick_right()
			thisaxis.yaxis.set_label_position("right")
			thisaxis.set_xlabel('frames')
			thisaxis.get_yaxis().set_major_locator(MaxNLocator(prune='both'))
			thisaxis.get_xaxis().set_major_locator(MaxNLocator(prune='lower'))
			if p == max(appor)+1:
				thisaxis.set_xlabel('$\mathsf{\sigma_a,\sigma_b(nm^{-1})}$',fontsize=14)
			else:
				thisaxis.set_xticklabels([])		
	plt.subplots_adjust(hspace = 0)
	plt.subplots_adjust(wspace = 0)
	plt.show()	
	
#---Advanced plotting method
if do_stacked_plot:
	fig = plt.figure(figsize=(10,10))
	gs = gridspec.GridSpec(max(appor)+1,3,wspace=0.0,hspace=0.0)
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
		ax.set_ylim(0,1.1*maxpeak)
		ax.axvline(x=0,ls='-',lw=1,c='k')
		ax.legend(loc=2,prop={'size':12})
		ax.set_ylabel('frequency')
		ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
		ax.grid(True)
		if a == len(axes_maxcurv)-1:
			ax.spines['bottom'].set_position(('outward', 10))
			ax.set_xlabel('$\mathsf{H_{max}(nm^{-1})}$',fontsize=16)
			second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
			ax.spines['second_bottom'] = second_bottom
			ax.set_xticks(arange(-0.10,0.12,0.02))
		else:
			ax.set_xticks(arange(-0.10,0.12,0.02))
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
			ax.set_title(r'$\textbf{extents of curvature}$')
		ax.grid(True)
		#ax.yaxis.tick_right()
		#ax.yaxis.set_label_position("right")
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
		ax.set_xticks(arange(5,30,5))
		if a == len(axes_sigmas)-1:
			ax.set_xlabel('$\mathsf{\sigma_a,\sigma_b(nm^{-1})}$',fontsize=14)
		else:
			ax.set_xticklabels([])		
	#---plot areas
	for p in range(len(analysis_descriptors[analysis_plan])):
		if appor[p] > 0 and appor[p] == appor[p-1]:
			thisaxis = axes_maxcurv[-1]
		else:
			thisaxis = fig.add_subplot(gs[appor[p],3])
		area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
		posarea = array([area_per_tile*area_counts[i,2] for i in range(len(area_counts))])
		negarea = array([area_per_tile*area_counts[i,3] for i in range(len(area_counts))])



	plt.subplots_adjust(hspace = 0)
	plt.subplots_adjust(wspace = 0)
	plt.show()	


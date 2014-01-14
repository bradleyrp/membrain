#!/usr/bin/python -i

if 1:

	from membrainrunner import *

	import os
	from scipy.optimize import leastsq
	import matplotlib as mpl
	from pylab import *

	skip = 1
	framecount = None
	location = ''
	execfile('locations.py')
	
	which_brewer_colors = [0,1,2,3,4,5,6,7]
	clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
	clrs2 = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
	mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	mpl.rcParams['axes.linewidth'] = 2.0

	
	#---analysis plan
	analysis_descriptors = [
		('pkl.dimple.v614-stress.md.part0002.rerun.pkl',),
		('pkl.dimple.v612-stress.md.part0003.pkl',),
		('pkl.dimple.v550.md.part0006.300000-400000-200.pkl',),
		('pkl.dimple.v550.md.part0006.300000-400000-200.testshift10.pkl',),
		('pkl.dimple.v550.md.part0006.300000-400000-200.testshift01.pkl',),
		('pkl.dimple.v550.md.part0006.300000-400000-200.testshift11.pkl',),
		('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',),
		('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',)]

	do_single_plot = True
	do_stacked_plot = True
	do_opposite_signs = False

	analysis_plan = range(len(analysis_descriptors))
	names = ('ENTHx4','ENTHx1','control A','control B','control C','control D',
		'exo70x2 (anti)','exo70x2 (parallel)')
	appor = (0,1,2,2,2,2,3,4)
	ccodes = [(clrs[0],clrs[1]),(clrs[2],clrs[3]),
		(clrs[4],clrs[5]),(clrs[6],clrs[7]),(clrs[0],clrs[1]),(clrs[2],clrs[3]),
		(clrs[0],clrs[1]),(clrs[2],clrs[3])]
	fillcodes = (1,1,0,0,0,0,1,1)
	
	
	'''
	possible plots
	1. all ENTH systems on one plot
	2. all ENTH systems on a plot stacked above Exo70 systems
	3. all systems on one plot
	'''
	
	results_stack = []
	for pnum in analysis_plan:
		results_stack.append(pickle.load(open(pickles+analysis_descriptors[pnum][0],'r')))
	
if do_single_plot:
	fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(10,4))
	nbins = 20
	minval,maxval = -0.10,0.10
	maxpeak = 0
	for p in range(len(analysis_plan)):
		pnum = analysis_plan[p]
		expected_direction = results_stack[pnum][0].notes[([i[0] 
			for i in results_stack[pnum][0].notes].index('expected_direction'))][1]
		order = ((1,0) if expected_direction == 1 else (0,1))
		for o in order:
			params = results_stack[pnum][o].get(['type','params'])
			maxhs = results_stack[pnum][o].get(['type','maxhs'])
			maxhxys = results_stack[pnum][o].get(['type','maxhxys'])
			validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
			if max(hist0) > maxpeak: maxpeak = max(hist0)
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			if o == 1:
				axes.plot(mid0,hist0,c=ccodes[p][o],alpha=1.,lw=2,label=names[p])
			else:
				axes.plot(mid0,hist0,c=ccodes[p][o],alpha=1.,lw=2)
			if fillcodes[p]:
				axes.fill_between(mid0,hist0,[0 for i in mid0],facecolor=ccodes[p][o],alpha=0.2,
					interpolate=True)
	axes.axvline(x=0,ls='-',lw=1,c='k')
	plt.show()	


if do_stacked_plot:
	fig, axes = plt.subplots(nrows=max(appor)+1,ncols=1,figsize=(6,2*(max(appor)+1)))
	nbins = 20
	minval,maxval = -0.10,0.10
	maxpeak = 0
	for p in range(len(analysis_plan)):
		pnum = analysis_plan[p]
		expected_direction = results_stack[pnum][0].notes[([i[0] 
			for i in results_stack[pnum][0].notes].index('expected_direction'))][1]
		order = ((1,0) if expected_direction == 1 else (0,1))
		for o in order:
			params = results_stack[pnum][o].get(['type','params'])
			maxhs = results_stack[pnum][o].get(['type','maxhs'])
			maxhxys = results_stack[pnum][o].get(['type','maxhxys'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			if o == 1:
				axes[appor[p]].plot(mid0,hist0,c=ccodes[p][o],alpha=1.,lw=2,label=names[p])
			else:
				axes[appor[p]].plot(mid0,hist0,c=ccodes[p][o],alpha=1.,lw=2)
			if fillcodes[p]:
				axes[appor[p]].fill_between(mid0,hist0,[0 for i in mid0],facecolor=ccodes[p][o],alpha=0.2,
					interpolate=True)
			if max(hist0) > maxpeak: maxpeak = max(hist0)
	for p in range(len(appor)):
		axes[appor[p]].set_ylim(0,1.1*maxpeak)
		axes[appor[p]].axvline(x=0,ls='-',lw=1,c='k')
		axes[appor[p]].legend()
	plt.show()	
		


if 0:
	means='linear'
	scales='linear'
	fullrange=True
	extrasuffix=''
	sigma_calc='mode'


	single_plot = 1
	if single_plot:

		params_neg = result_data_collection[0].get(['type','params'])
		maxhs_neg = result_data_collection[0].get(['type','maxhs'])
		maxhxys_neg = result_data_collection[0].get(['type','maxhxys'])

		params_poz = result_data_collection[1].get(['type','params'])
		maxhs_poz = result_data_collection[1].get(['type','maxhs'])
		maxhxys_poz = result_data_collection[1].get(['type','maxhxys'])
	
		params_batch = (params_neg,params_poz)
		maxhs_batch = (maxhs_neg,maxhs_poz)
		maxhxys_batch = (maxhxys_neg,maxhxys_poz)

		'''
		params = (params_neg,params_poz)
		maxhs = (maxhs_neg,maxhs_poz)
		maxhxys = (maxhxys_neg,maxhxys_poz)

		(params_neg,params_poz) = params
		(maxhs_neg,maxhs_poz) = maxhs
		(maxhxys_neg,maxhxys_poz) = maxhxys
		'''
	
		#---top plot
		fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(10,4))
		nbins = 20
		minval,maxval = -0.10,0.10
		
		params, maxhs, maxhxys = params_neg, maxhs_neg, maxhxys_neg
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
		#---nanometer correction
		validhs = [10*maxhs[i] for i in validhis]
		hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		axes.plot(mid0,hist0,c=clrs[1],alpha=1.,lw=2)
		axes.fill_between(mid0,hist0,[0 for i in mid0],facecolor=clrs[1],alpha=0.35,interpolate=True)
		
		params, maxhs, maxhxys = params_poz, maxhs_poz, maxhxys_poz
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
		#---nanometer correction
		validhs = [10*maxhs[i] for i in validhis]
		#---plot
		hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		axes.plot(mid0,hist0,c=clrs[0],alpha=1.,lw=2)
		axes.fill_between(mid0,hist0,[0 for i in mid0],facecolor=clrs[0],alpha=0.35,interpolate=True)
	


		plt.show()



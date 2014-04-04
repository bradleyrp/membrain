#!/usr/bin/python

'''
This script is the classic version of the second-step (reproc) of the dimple-fitting procedure, raised from 
the dead in order to rectify the old and new versions. It was formerly known as 
script-reproc-dimple-backup-2014.03.06.py script-reproc-dimple-classic.py
'''

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
		r'$\textbf{{ENTH}\ensuremath{\times}4}$',1,None),
	('pkl.dimple.v612-stress.md.part0003.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{ENTH}\ensuremath{\times}1}$',1,None),
	('pkl.dimple.v612-stress.md.part0003.prot-v614.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{ENTH}\ensuremath{\times}1}$',1,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.invert.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control, invert}$',0,'inverttest'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v614.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v700.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v700.pkl',(clrs[4],clrs[5]),
		r'$\textbf{control shift}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v700.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control}$',0,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-16-16.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+16,16}$',0,(16,16)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-8-8.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+8,8}$',0,(8,8)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-24-24.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+24,24}$',0,(24,24)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-0-16.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+0,16}$',0,(0,16)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-16-0.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+16,0}$',0,(16,0)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-2-2.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+2,2}$',0,(2,2)),
    ('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-0-2.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+0,2}$',0,(0,2)),
    ('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-2-0.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+2,0}$',0,(2,0)),
    ('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.shift-4-4.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+4,4}$',0,(4,4)),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.phasetest.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control+rand}$',0,'phasetest'),
	('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (para)}}$',1,None),
	('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',(clrs[2],clrs[3]),
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (anti)}}$',1,None),
	('pkl.dimple.v550.md.part0006.300000-400000-200.dummytest.pkl',(clrs[6],clrs[7]),
		r'$\textbf{control, invert}$',0,None),
	('pkl.dimple.v614.s9-lonestar.md.part0004.120000-220000-200.s9.120000-220000-200.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4\,s9}$',1,None),
	('pkl.dimple.v614.s6-sim-lonestar.md.part0002.40000-140000-200.s6.40000-140000-200.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4\,s6}$',1,None),
	('pkl.dimple.v614-stress.md.part0002.rerun.stress.part0002.rerun.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4\,p2}$',1,None)]

#---methods
do_stacked_plot = True
do_hmax_vs_sigmas = False
do_errorlook = False
do_errorlook_sigma = False
do_sigma_vs_hmax = False
do_resid_1d = False
do_pubplot = False
		
#---pre-made analysis routines
plotspecs = 'enthpub'
if plotspecs == 'both':
	analysis_plan = slice(None,-1)
	figoutnamebase = 'fig-dimple-master-summary-ENTH-EXO70-classic-check2'
	figsize = (14,16)
	analysis_plan = [0,1,2,3,4,9,10,11,12,13,14,15,16,17,18,19,20,23]
	analysis_plan = [0,1,2,3,4]
	analysis_plan = [0,1,2,3,4,9,10,11,12,13,23,24,25]
	analyses = [analysis_descriptors[i] for i in analysis_plan]
	appor = range(len(analysis_plan))
elif plotspecs == 'enth':
	analysis_plan = slice(0,3)
	appor = (0,1,2)
	figoutnamebase = 'fig-dimple-master-summary-ENTH'
	figsize = (14,8)
elif plotspecs == 'enthpub':
	analysis_plan = [0,2,3]
	analyses = [analysis_descriptors[i] for i in analysis_plan]
	do_stacked_plot = False
	do_pubplot = True
	figsize = (10,4)
	figoutnamebase = 'fig-dimple-manuscript-ENTH'

#---method settings
do_resid_filter = True
resid_filter = 6.
show_areas_on_master_plot = True

#---parameters
nbins = 20
nbins_sigma = 20
minval,maxval = -0.10,0.10
minval_sigma,maxval_sigma = 0,30
maxhfilter = [0.001,0.1]
maxhrange = (-0.06,0.06)
maxhstep = 0.01
sigmamax = 30
which_orders = [2]
which_orders_extents = 2

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	z0 = 0
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)
		
extra_locs = [
	'./','structures-broken-transposer-error/',
	'backup-2013.20.21-enth-review-pickles/']
def unpickle_special(name):
	'''Custom un-pickle code to search for previous pickles.'''
	mset = None
	for loc in extra_locs:
		mset = unpickle(pickles+loc+name)
		if mset != None: break
	print mset
	return mset

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
results_stack = []
for pnum in range(len(analyses)):
	results_stack.append(pickle.load(open(pickles+analyses[pnum][0],'r')))
	#---hack to fix plot on v614 additions
	if analyses[pnum][0] in [
		'pkl.dimple.v614.s9-lonestar.md.part0004.120000-220000-200.s9.120000-220000-200.pkl',
		'pkl.dimple.v614.s6-sim-lonestar.md.part0002.40000-140000-200.s6.40000-140000-200.pkl',
		'pkl.dimple.v614-stress.md.part0002.rerun.stress.part0002.rerun.pkl']:
		for j in range(3):
			results_stack[-1][j].notes[([i[0] 
				for i in results_stack[-1][0].notes].index('expected_direction'))][1] = 1

#---publication-style plots, ENTH results
if do_pubplot:
	#---this is custom plot for IET manuscript
	#---note overall grid is really 3 rows (systems) and 5 columns (3 for curvs, 2 for area)
	#---plots
	fig = plt.figure(figsize=figsize)
	ax_hmax = []
	ax_sigs = []
	ax_area = []
	ax_area_distn = []
	#---upper left has ENTH systems with Hmax and sigma
	gs00 = gridspec.GridSpec(2,3)
	gs00.update(left=0.02, right=0.58,top=0.98,bottom=0.333+0.02,wspace=0.0,hspace=0.0)
	ax = plt.subplot(gs00[0,0:2]);ax_hmax.append(ax)
	ax = plt.subplot(gs00[1,0:2]);ax_hmax.append(ax)
	ax = plt.subplot(gs00[0,2]);ax_sigs.append(ax)
	ax = plt.subplot(gs00[1,2]);ax_sigs.append(ax)
	#---lower left has control systems with Hmax and sigma
	gs10 = gridspec.GridSpec(1,3)
	gs10.update(left=0.02, right=0.58,top=0.333-0.02,bottom=0.02,wspace=0.0,hspace=0.0)
	ax = plt.subplot(gs10[0,0:2]);ax_hmax.append(ax)
	ax = plt.subplot(gs10[0,2]);ax_sigs.append(ax)
	#---upper right has ENTH systems with areas
	gs01 = gridspec.GridSpec(2,2)
	gs01.update(left=0.62,right=0.98,top=0.98,bottom=0.333+0.06,wspace=0.0,hspace=0.0)
	ax = plt.subplot(gs01[0,0]);ax_area.append(ax)
	ax = plt.subplot(gs01[0,1]);ax_area_distn.append(ax)
	ax = plt.subplot(gs01[1,0]);ax_area.append(ax)
	ax = plt.subplot(gs01[1,1]);ax_area_distn.append(ax)
	#---lower right has control systems with areas
	gs11 = gridspec.GridSpec(1,2)
	gs11.update(left=0.62,right=0.98,top=0.333-0.06,bottom=0.02,wspace=0.0,hspace=0.0)	
	ax = plt.subplot(gs11[0,0]);ax_area.append(ax)
	ax = plt.subplot(gs11[0,1]);ax_area_distn.append(ax)
	#---plot maximum mean curvatures	
	maxpeak = 0
	for p in range(len(analyses)):
		print 'p = '+str(p)
		thisaxis = ax_hmax[p]
		ccodes = analyses[p][1]
		name = analyses[p][2]
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		meanslist = [0,0,0]
		valid_frames = [0,0,0]
		for o in which_orders:
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])

			if do_resid_filter:
				target_zones = results_stack[p][order[o]].get(['type','target_zones'])
				resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
					for j in range(len(maxhs))]
				validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
					and abs(10*maxhs[i]) < maxhfilter[1]) and resids[i] < resid_filter]
			else:			
				validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
					and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			print 'mean = '+str(mean(validhs))
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=False,weights=[1./len(validhs) 
				for i in validhs],range=(minval,maxval))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			if o == 1:
				#thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2,label=name)
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2)
			elif o == 0:
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2)
			else:
				thisaxis.plot(mid0,hist0,'-',c=ccodes[1],alpha=1.,lw=2)
			thisaxis.fill_between(mid0,hist0,[0 for i in mid0],facecolor=ccodes[1],
				alpha=0.2,interpolate=True)
			if max(hist0) > maxpeak: maxpeak = max(hist0)
			meanslist[order[o]] = mean(validhs)
			valid_frames[order[o]] = len(validhis)
		if 0:
			textline = r'$\left\langle H_{max}\right\rangle =\textrm{('+str('%3.3f'%meanslist[0])+','+\
				str('%3.3f'%meanslist[1])+','+str('%3.3f'%meanslist[2])+')}$'+\
				'\n['+str('%d'%valid_frames[0])+','+str('%d'%valid_frames[1])+\
				','+str('%d'%valid_frames[2])+'] of '+str(len(maxhs))+'  '
			thisaxis.text(0.98,0.9,textline,transform=thisaxis.transAxes,fontsize=12,
				horizontalalignment='right',verticalalignment='top')
		thisaxis.set_ylabel(name,rotation='horizontal')
	for a in range(len(ax_hmax)):
		ax = ax_hmax[a]
		if a == 0:
			ax.set_title(r'$\textbf{curvature}$',fontsize=16)
		ax.set_ylim(0,1.2*maxpeak)
		ax.axvline(x=0,ls='-',lw=1,c='k')
		ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
		ax.grid(True)
		ax.set_xlim(maxhrange)
		ax.set_yticklabels([])
		if a == len(ax_hmax)-1:
			ax.set_xticks(arange(maxhrange[0],maxhrange[1]+0.001,maxhstep))
			ax.spines['bottom'].set_position(('outward', 10))
			ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
			second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
			ax.spines['second_bottom'] = second_bottom
		else:
			ax.set_xticks(arange(maxhrange[0],maxhrange[1]+.0001,maxhstep))
			ax.set_xticklabels([])
	#---plot extents of curvature
	axes_sigmas = []
	maxpeak = 0
	sigmeans = []
	sigmodes = []
	valid_frames = []
	total_frames = []
	for p in range(len(analyses)):
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		if which_orders_extents != None:
			order = [which_orders_extents]
		else:
			order = [1] if expected_direction == 1 else [0]
		o = 0
		thisaxis = ax_sigs[p]
		params = results_stack[p][order[o]].get(['type','params'])
		maxhs = results_stack[p][order[o]].get(['type','maxhs'])
		ccodes = analyses[p][1]
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
			and abs(10*maxhs[i]) < maxhfilter[1])]
		sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][4])/10. < sigmamax]
		sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][5])/10. < sigmamax]
		hist0,binedge0 = numpy.histogram(sigma_x,bins=nbins_sigma,normed=True,density=True,
			range=(minval_sigma,maxval_sigma))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		thisaxis.plot(mid0,hist0,c=ccodes[1],alpha=1.,lw=2)
		hist1,binedge1 = numpy.histogram(sigma_y,bins=nbins_sigma,normed=True,density=True,
			range=(minval_sigma,maxval_sigma))
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		thisaxis.plot(mid1,hist1,c=ccodes[0],alpha=1.,lw=2)
		if max(max(hist0),max(hist1)) > maxpeak: maxpeak = max(max(hist0),max(hist1))
		sigmeans.append([mean(sigma_x),mean(sigma_y)])
		sigmodes.append([mid0[argmax(hist0)],mid1[argmax(hist1)]])
		valid_frames.append([len(sigma_x),len(sigma_y)])
		total_frames.append(len(maxhs))
	for a in range(len(ax_sigs)):
		ax = ax_sigs[a]
		if a == 0:
			ax.set_title(r'$\textbf{extent}$',fontsize=16)
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
		ax.set_xticks(arange(5,maxval_sigma+0.001,5))
		if a == len(ax_sigs)-1:
			ax.set_xlabel('$\mathsf{\sigma_a,\sigma_b\,(nm)}$',fontsize=14)
		else:
			ax.set_xticklabels([])		
		if 0:
			textline = 'means: ['+str('%3.1f'%sigmeans[a][0])+','+str('%3.1f'%sigmeans[a][1])+\
				']\n'+'modes: ['+str('%3.1f'%sigmodes[a][0])+','+str('%3.1f'%sigmodes[a][1])+']'+\
				'\n['+str('%d'%valid_frames[a][0])+','+str('%d'%valid_frames[a][1])+'] of '+\
				str(total_frames[a])+'  '
			ax.text(0.98,0.9,textline,transform=ax.transAxes,fontsize=12,horizontalalignment='right',
				verticalalignment='top')	
	#---plot areas
	maxpeak = 0
	numframes = []
	for p in range(len(analyses)):
		thisaxis = ax_area[p]
		mset = unpickle_special(results_stack[p][2].getnote('startpickle'))
		vecs = mean(mset.vecs,axis=0)
		area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
		target_zones = results_stack[p][2].get(['type','target_zones'])
		surfpbc = []
		if type(analyses[p][4]) == str and analyses[p][4] == 'phasetest':
			for fr in range(len(target_zones)):
				randshift = [randint(0,31) for i in range(2)]
				surfpbc.append([[mset.surf[fr][(i+randshift[0])%(mset.griddims[0]-2)]\
					[(j+randshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))])
		elif type(analyses[p][4]) == str and analyses[p][4] == 'inverttest':
			surfpbc = -1*array(mset.surf)
		elif type(analyses[p][4]) == tuple:
			fixedshift = analyses[p][4]
			for fr in range(len(target_zones)):
				surfpbc.append([[mset.surf[fr][(i+fixedshift[0])%(mset.griddims[0]-2)]\
				[(j+fixedshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
				for i in range(1*(mset.griddims[0]-1))])
		else:
			surfpbc = mset.surf
		#---hack to fix weird indexing not sure why this script isn't executing perfectly
		nareas = min([len(target_zones),len(surfpbc)])
		poz_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
			[int(i[1]*(mset.griddims[1]-1)/vecs[1])]>0. for i in target_zones[j]]) 
			for j in range(nareas)]
		neg_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
			[int(i[1]*(mset.griddims[1]-1)/vecs[1])]<0. for i in target_zones[j]]) 
			for j in range(nareas)]
		posarea = array([area_per_tile*poz_area_counts[i] for i in range(len(poz_area_counts))])
		negarea = array([area_per_tile*neg_area_counts[i] for i in range(len(neg_area_counts))])
		thisaxis.plot(posarea,'r-',label='$z>0$',lw=1,alpha=0.8)
		thisaxis.plot(negarea,'b-',label='$z<0$',lw=1,alpha=0.8)
		t = range(len(poz_area_counts))
		thisaxis.fill_between(t, negarea,posarea, facecolor='b',alpha=0.25,where=negarea>posarea)
		thisaxis.fill_between(t, posarea,negarea, facecolor='r',alpha=0.25,where=negarea<posarea)
		if max(max(posarea),max(negarea)) > maxpeak: maxpeak = max(max(posarea),max(negarea))
		numframes.append(len(target_zones))
	for a in range(len(ax_area)):
		ax = ax_area[a]
		if a == 0 and False:
			ax.set_title(r'$\textbf{areas (+/-)}$')
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
		ax.set_xlim((0,numframes[a]))		
		if a == len(ax_area)-1:
			ax.set_xlabel('frame',fontsize=14)
		#---modification for extra axis here
		elif a != 1:
			ax.set_xticklabels([])		
	#---plot area histograms
	minval_areas = 0
	maxval_areas = maxpeak
	maxfreq = 0
	for p in range(len(analyses)):
		thisaxis = ax_area_distn[p]
		mset = unpickle_special(results_stack[p][2].getnote('startpickle'))
		vecs = mean(mset.vecs,axis=0)
		area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
		target_zones = results_stack[p][2].get(['type','target_zones'])
		surfpbc = []
		if type(analyses[p]) == str and analyses[p][4] == 'phasetest':
			for fr in range(len(target_zones)):
				randshift = [randint(0,31) for i in range(2)]
				surfpbc.append([[mset.surf[fr][(i+randshift[0])%(mset.griddims[0]-2)]\
					[(j+randshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))])
		elif type(analyses[p][4]) == str and analyses[p][4] == 'inverttest':
			surfpbc = -1*array(mset.surf)
		elif type(analyses[p][4]) == tuple:
			fixedshift = analyses[p][4]
			for fr in range(len(target_zones)):
				surfpbc.append([[mset.surf[fr][(i+fixedshift[0])%(mset.griddims[0]-2)]\
				[(j+fixedshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
				for i in range(1*(mset.griddims[0]-1))])
		else:
			surfpbc = mset.surf
		poz_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
			[int(i[1]*(mset.griddims[1]-1)/vecs[1])]>0. for i in target_zones[j]]) 
			for j in range(len(target_zones))]
		neg_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
			[int(i[1]*(mset.griddims[1]-1)/vecs[1])]<0. for i in target_zones[j]]) 
			for j in range(len(target_zones))]
		posarea = array([area_per_tile*poz_area_counts[i] for i in range(len(poz_area_counts))])
		negarea = array([area_per_tile*neg_area_counts[i] for i in range(len(neg_area_counts))])
		hist0,binedge0 = numpy.histogram(posarea,bins=nbins_sigma,normed=False,
			weights=[1./len(negarea) for i in negarea],
			range=(minval_areas,maxval_areas))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		thisaxis.plot(hist0,mid0,c='r',alpha=1.,lw=1,label='$z>0$')
		hist1,binedge1 = numpy.histogram(negarea,bins=nbins_sigma,normed=False,
			weights=[1./len(negarea) for i in negarea],
			range=(minval_areas,maxval_areas))
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		thisaxis.plot(hist1,mid1,c='b',alpha=1.,lw=1,label='$z<0$')
		thisaxis.fill_betweenx(mid0,[0 for i in mid0],hist0,facecolor='r',
			alpha=0.2)
		thisaxis.fill_betweenx(mid1,[0 for i in mid1],hist1,facecolor='b',
			alpha=0.2)
		if max(max(hist0),max(hist1)) > maxfreq: maxfreq = max(max(hist0),max(hist1))
	for a in range(len(ax_area_distn)):
		ax = ax_area_distn[a]
		#---disabled
		if a == 0 and False:
			ax.set_title(r'$\textbf{areas (+/-)}$')
		ax.grid(True)
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		ax.set_ylabel(r'$\textbf{area}$'+'\n'+r'$\textbf{(nm\ensuremath{{}^{2}})}$',fontsize=10)
		maxfreq = 0.5
		ax.set_xlim(0,1.0*maxfreq)
		ax.get_yaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
		#---disabled all frequency labels
		if a == len(ax_area_distn)-1 and False:
			ax.set_xlabel('frequency',fontsize=14)
		else:
			ax.set_xticklabels([])
	#---extra frame counter since frame counts different
	ax = ax_area[1]
	ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
	ax.set_xlabel('frame')
	#---shared area label
	fig.text(0.78,0.995,r'$\textbf{fitted areas}$',ha='center',va='bottom',rotation='horizontal',fontsize=16)
	#---extra legend for areas
	outsidelegendax = ax_area_distn[-1]
	axpos = outsidelegendax.get_position()
	outlegend = outsidelegendax.legend(loc='right',prop={'size':12,'weight':'bold'})
	for legobj in outlegend.legendHandles:
		legobj.set_linewidth(3.0)
	#---finalize	
	if do_resid_filter:
		figoutname = figoutnamebase+'-rmsd-filter-'+str(resid_filter)+'A-maxhfilter-'+\
			str(maxhfilter[0])+'-'+str(maxhfilter[1])+'.png'
	else:
		figoutname = figoutnamebase+'-maxhfilter-'+str(maxhfilter[0])+'-'+str(maxhfilter[1])+'.png'
	plt.savefig(pickles+figoutname,dpi=500,bbox_inches='tight')
	plt.show()
	
#---summary plots
if do_stacked_plot:
	fig = plt.figure(figsize=figsize)
	gs = gridspec.GridSpec(max(appor)+1,(6 if show_areas_on_master_plot else 4),wspace=0.0,hspace=0.0)
	#---plot maximum mean curvatures	
	maxpeak = 0
	axes_maxcurv = []
	for p in range(len(analyses)):
		if appor[p] > 0 and appor[p] == appor[p-1]:
			thisaxis = axes_maxcurv[-1]
		else:
			thisaxis = fig.add_subplot(gs[appor[p],0:2])
		ccodes = analyses[p][1]
		name = analyses[p][2]
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		meanslist = [0,0,0]
		valid_frames = [0,0,0]
		for o in range(len(order)):
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			if do_resid_filter:
				target_zones = results_stack[p][order[o]].get(['type','target_zones'])
				resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
					for j in range(len(maxhs))]
				validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
					and abs(10*maxhs[i]) < maxhfilter[1]) and resids[i] < resid_filter]
			else:			
				validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
					and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=False,weights=[1./len(validhs) 
				for i in validhs],range=(minval,maxval))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			if o == 1:
				#thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2,label=name)
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2)
			elif o == 0:
				thisaxis.plot(mid0,hist0,'-',c=ccodes[o],alpha=1,lw=2)
			else:
				thisaxis.plot(mid0,hist0,'--',c='k',alpha=1.,lw=2)
			if fillcode and o != 2:
				thisaxis.fill_between(mid0,hist0,[0 for i in mid0],facecolor=ccodes[o],
					alpha=0.2,interpolate=True)
			if max(hist0) > maxpeak: maxpeak = max(hist0)
			meanslist[order[o]] = mean(validhs)
			valid_frames[order[o]] = len(validhis)
		textline = r'$\left\langle H_{max}\right\rangle =\textrm{('+str('%3.3f'%meanslist[0])+','+\
			str('%3.3f'%meanslist[1])+','+str('%3.3f'%meanslist[2])+')}$'+\
			'\n['+str('%d'%valid_frames[0])+','+str('%d'%valid_frames[1])+\
			','+str('%d'%valid_frames[2])+'] of '+str(len(maxhs))+'  '
		thisaxis.text(0.98,0.9,textline,transform=thisaxis.transAxes,fontsize=12,
			horizontalalignment='right',verticalalignment='top')
		if not (appor[p] > 0 and appor[p] == appor[p-1]):
			axes_maxcurv.append(thisaxis)
		thisaxis.set_ylabel(name,rotation='horizontal')
	for a in range(len(axes_maxcurv)):
		ax = axes_maxcurv[a]
		if a == 0:
			ax.set_title(r'$\textbf{mean curvatures}$')
		ax.set_ylim(0,1.2*maxpeak)
		ax.axvline(x=0,ls='-',lw=1,c='k')
		ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
		ax.grid(True)
		ax.set_xlim(maxhrange)
		ax.set_yticklabels([])
		if a == len(axes_maxcurv)-1:
			ax.set_xticks(arange(maxhrange[0],maxhrange[1]+0.001,maxhstep))
			ax.spines['bottom'].set_position(('outward', 10))
			ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
			second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
			ax.spines['second_bottom'] = second_bottom
		else:
			ax.set_xticks(arange(maxhrange[0],maxhrange[1]+.0001,maxhstep))
			ax.set_xticklabels([])
	#---plot extents of curvature
	axes_sigmas = []
	maxpeak = 0
	sigmeans = []
	sigmodes = []
	valid_frames = []
	total_frames = []
	for p in range(len(analyses)):
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = [1] if expected_direction == 1 else [0]
		o = 0
		thisaxis = plt.subplot(gs[appor[p],2])
		params = results_stack[p][order[o]].get(['type','params'])
		maxhs = results_stack[p][order[o]].get(['type','maxhs'])
		ccodes = analyses[p][1]
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
			and abs(10*maxhs[i]) < maxhfilter[1])]
		sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][4])/10. < sigmamax]
		sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][5])/10. < sigmamax]
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
		sigmeans.append([mean(sigma_x),mean(sigma_y)])
		sigmodes.append([mid0[argmax(hist0)],mid1[argmax(hist1)]])
		valid_frames.append([len(sigma_x),len(sigma_y)])
		total_frames.append(len(maxhs))
	for a in range(len(axes_sigmas)):
		ax = axes_sigmas[a]
		if a == 0:
			ax.set_title(r'$\textbf{extents (+/-)}$')
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
		ax.set_xticks(arange(5,maxval_sigma+0.001,5))
		if a == len(axes_sigmas)-1:
			ax.set_xlabel('$\mathsf{\sigma_a,\sigma_b\,(nm)}$',fontsize=14)
		else:
			ax.set_xticklabels([])		
		textline = 'means: ['+str('%3.1f'%sigmeans[a][0])+','+str('%3.1f'%sigmeans[a][1])+\
			']\n'+'modes: ['+str('%3.1f'%sigmodes[a][0])+','+str('%3.1f'%sigmodes[a][1])+']'+\
			'\n['+str('%d'%valid_frames[a][0])+','+str('%d'%valid_frames[a][1])+'] of '+\
			str(total_frames[a])+'  '
		ax.text(0.98,0.9,textline,transform=ax.transAxes,fontsize=12,horizontalalignment='right',
			verticalalignment='top')
	#---plot extents of curvature, for the basic shadow
	axes_sigmas = []
	maxpeak = 0
	sigmeans = []
	sigmodes = []
	valid_frames = []
	total_frames = []
	for p in range(len(analyses)):
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = [1] if expected_direction == 1 else [0]
		order = [2]
		o = 0
		thisaxis = plt.subplot(gs[appor[p],3])
		params = results_stack[p][order[o]].get(['type','params'])
		maxhs = results_stack[p][order[o]].get(['type','maxhs'])
		ccodes = analyses[p][1]
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > maxhfilter[0] 
			and abs(10*maxhs[i]) < maxhfilter[1])]
		sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][4])/10. < sigmamax]
		sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
			and abs(params[i][5])/10. < sigmamax]
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
		sigmeans.append([mean(sigma_x),mean(sigma_y)])
		sigmodes.append([mid0[argmax(hist0)],mid1[argmax(hist1)]])
		valid_frames.append([len(sigma_x),len(sigma_y)])
		total_frames.append(len(maxhs))
	for a in range(len(axes_sigmas)):
		ax = axes_sigmas[a]
		if a == 0:
			ax.set_title(r'$\textbf{extents}$')
		ax.grid(True)
		ax.set_ylim(0,1.1*maxpeak)
		ax.set_yticklabels([])		
		ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
		ax.set_xticks(arange(5,maxval_sigma+0.001,5))
		if a == len(axes_sigmas)-1:
			ax.set_xlabel('$\mathsf{\sigma_a,\sigma_b\,(nm)}$',fontsize=14)
		else:
			ax.set_xticklabels([])		
		textline = 'means: ['+str('%3.1f'%sigmeans[a][0])+','+str('%3.1f'%sigmeans[a][1])+\
			']\n'+'modes: ['+str('%3.1f'%sigmodes[a][0])+','+str('%3.1f'%sigmodes[a][1])+']'+\
			'\n['+str('%d'%valid_frames[a][0])+','+str('%d'%valid_frames[a][1])+'] of '+\
			str(total_frames[a])+'  '
		ax.text(0.95,0.9,textline,transform=ax.transAxes,fontsize=12,horizontalalignment='right',
			verticalalignment='top')

	if show_areas_on_master_plot:
		#---plot areas
		axes_areas = []
		maxpeak = 0
		for p in range(len(analyses)):
			if appor[p] > 0 and appor[p] == appor[p-1]:
				thisaxis = axes_areas[-1]
				repeat = True
			else:
				thisaxis = fig.add_subplot(gs[appor[p],4])
				repeat = False
			mset = unpickle_special(results_stack[p][2].getnote('startpickle'))
			vecs = mean(mset.vecs,axis=0)
			area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
			target_zones = results_stack[p][2].get(['type','target_zones'])
			surfpbc = []
			if type(analyses[p][4]) == str and analyses[p][4] == 'phasetest':
				for fr in range(len(target_zones)):
					randshift = [randint(0,31) for i in range(2)]
					surfpbc.append([[mset.surf[fr][(i+randshift[0])%(mset.griddims[0]-2)]\
						[(j+randshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
						for i in range(1*(mset.griddims[0]-1))])
			elif type(analyses[p][4]) == str and analyses[p][4] == 'inverttest':
				surfpbc = -1*array(mset.surf)
			elif type(analyses[p][4]) == tuple:
				fixedshift = analyses[p][4]
				for fr in range(len(target_zones)):
					surfpbc.append([[mset.surf[fr][(i+fixedshift[0])%(mset.griddims[0]-2)]\
					[(j+fixedshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))])
			else:
				surfpbc = mset.surf
			poz_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
				[int(i[1]*(mset.griddims[1]-1)/vecs[1])]>0. for i in target_zones[j]]) 
				for j in range(len(target_zones))]
			neg_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
				[int(i[1]*(mset.griddims[1]-1)/vecs[1])]<0. for i in target_zones[j]]) 
				for j in range(len(target_zones))]
			posarea = array([area_per_tile*poz_area_counts[i] for i in range(len(poz_area_counts))])
			negarea = array([area_per_tile*neg_area_counts[i] for i in range(len(neg_area_counts))])
			thisaxis.plot(posarea,'r-',label=(None if repeat else '$z>0$' ),lw=1,alpha=0.8)
			thisaxis.plot(negarea,'b-',label=(None if repeat else '$z<0$'),lw=1,alpha=0.8)
			t = range(len(poz_area_counts))
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
				axpos = ax.get_position()
				outsidelegendax = ax
			else:
				ax.set_xticklabels([])		
		#---plot area histograms
		axes_area_hists = []
		minval_areas = 0
		maxval_areas = maxpeak
		maxfreq = 0
		for p in range(len(analyses)):
			if appor[p] > 0 and appor[p] == appor[p-1]:
				thisaxis = axes_area_hists[-1]
			else:
				thisaxis = fig.add_subplot(gs[appor[p],5])
			mset = unpickle_special(results_stack[p][2].getnote('startpickle'))
			vecs = mean(mset.vecs,axis=0)
			area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
			target_zones = results_stack[p][2].get(['type','target_zones'])
			surfpbc = []
			if type(analyses[p]) == str and analyses[p][4] == 'phasetest':
				for fr in range(len(target_zones)):
					randshift = [randint(0,31) for i in range(2)]
					surfpbc.append([[mset.surf[fr][(i+randshift[0])%(mset.griddims[0]-2)]\
						[(j+randshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
						for i in range(1*(mset.griddims[0]-1))])
			elif type(analyses[p][4]) == str and analyses[p][4] == 'inverttest':
				surfpbc = -1*array(mset.surf)
			elif type(analyses[p][4]) == tuple:
				fixedshift = analyses[p][4]
				for fr in range(len(target_zones)):
					surfpbc.append([[mset.surf[fr][(i+fixedshift[0])%(mset.griddims[0]-2)]\
					[(j+fixedshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))])
			else:
				surfpbc = mset.surf
			poz_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
				[int(i[1]*(mset.griddims[1]-1)/vecs[1])]>0. for i in target_zones[j]]) 
				for j in range(len(target_zones))]
			neg_area_counts = [sum([surfpbc[j][int(i[0]*(mset.griddims[0]-1)/vecs[0])]\
				[int(i[1]*(mset.griddims[1]-1)/vecs[1])]<0. for i in target_zones[j]]) 
				for j in range(len(target_zones))]
			posarea = array([area_per_tile*poz_area_counts[i] for i in range(len(poz_area_counts))])
			negarea = array([area_per_tile*neg_area_counts[i] for i in range(len(neg_area_counts))])
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
			ax.set_ylabel(r'$\textbf{area}$'+'\n'+r'$\textbf{(nm\ensuremath{{}^{2}})}$',fontsize=10)
			maxfreq = 0.5
			ax.set_xlim(0,1.0*maxfreq)
			ax.get_yaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
			ax.get_xaxis().set_major_locator(MaxNLocator(prune='both',nbins=4))
			if a == len(axes_area_hists)-1:
				ax.set_xlabel('frequency',fontsize=14)
			else:
				ax.set_xticklabels([])
		plt.subplots_adjust(hspace = 0)
		plt.subplots_adjust(wspace = 0)
	plt.tight_layout()
	axpos = outsidelegendax.get_position()
	outlegend = outsidelegendax.legend(loc='lower center',
		bbox_to_anchor=((axpos.x0+axpos.x1)/2.,axpos.y0/2.*0),
		ncol=2,prop={'size':12,'weight':'bold'},bbox_transform=gcf().transFigure)
	for legobj in outlegend.legendHandles:
		legobj.set_linewidth(3.0)
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
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
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][4])/10. < sigmamax]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][5])/10. < sigmamax]
			print 'mean hmax = '+str(mean(validhs))
			print 'mean sigmax = '+str(mean(sigma_x))
			print 'mean sigmay = '+str(mean(sigma_y))
			print 'no. frames = '+str(len(validhs))
	if do_resid_filter:
		figoutname = figoutnamebase+'-rmsd-filter-'+str(resid_filter)+'A-maxhfilter-'+\
			str(maxhfilter[0])+'-'+str(maxhfilter[1])+'.png'
	else:
		figoutname = figoutnamebase+'-maxhfilter-'+str(maxhfilter[0])+'-'+str(maxhfilter[1])+'.png'
	plt.savefig(pickles+figoutname,dpi=500,bbox_inches='tight')
	plt.show()	

#---compare the hmax and extent distributions
if do_hmax_vs_sigmas:
	ticknums = 5
	rounderx = 2
	roundery = 0
	fig = plt.figure(figsize=(12,4))
	gs = gridspec.GridSpec(3,len(analyses))
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		for o in order:
			ax = plt.subplot(gs[o,p])
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][4])/10. < sigmamax]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][5])/10. < sigmamax]
			meansig = [1./2*(abs(params[i][4])/10.+abs(params[i][5])/10.) 
				for i in validhis if len(shape(params[i])) > 0]
			H, xedges, yedges = histogram2d(validhs,meansig,bins=21,range=(maxhrange,(0,100)),
				weights=[1./150 for i in validhs])
			midx = (xedges[1:]+xedges[:-1])/2
			midy = (yedges[1:]+yedges[:-1])/2
			extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
			cmap = mpl.cm.jet
			cmap.set_bad(cmap(0),1.)
			ax.imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
				norm=None,cmap=cmap)
			ax.set_title(name,fontsize=7)
			ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=7)
			ax.set_ylabel('$\mathsf{\sigma_a,\sigma_b\,(nm)}$',fontsize=7)
			xts = [round(i,rounderx) for i in midx][::int(float(len(midx))/ticknums)]
			ax.axes.set_xticks([[round(i,rounderx) for i in midx].index(round(j,rounderx)) for j in xts])
			ax.axes.set_xticklabels(xts,fontsize=7)
			yts = [int(round(i,roundery)) for i in midy][::int(float(len(midx))/ticknums)]
			ax.axes.set_yticks([[round(i,roundery) for i in midy].index(round(j,roundery)) for j in yts])
			ax.axes.set_yticklabels(yts,fontsize=7)
	plt.savefig(pickles+'fig-dimple-hmax-vs-sigmas-filter-0.01-0.1.png',dpi=500,bbox_inches='tight')
	plt.show()

#---compare the hmax and residuals (2D histogram plots)
if do_errorlook:
	ticknums = 5
	rounderx = 2
	roundery = 0
	fig = plt.figure(figsize=(18,6))
	gs = gridspec.GridSpec(3,len(analyses),wspace=0.0,hspace=0.0)
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		order = (0,1,2)
		for o in order:
			ax = plt.subplot(gs[o,p])
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			target_zones = results_stack[p][order[o]].get(['type','target_zones'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
				for j in validhis]
			ax = plt.subplot(gs[o,p])
			H, xedges, yedges = histogram2d(validhs,resids,bins=21,range=(maxhrange,(0,20)),
				weights=[1./len(validhs) for i in validhs])
			midx = (xedges[1:]+xedges[:-1])/2
			midy = (yedges[1:]+yedges[:-1])/2
			extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
			cmap = mpl.cm.jet
			cmap.set_bad(cmap(0),1.)
			ax.imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
				norm=None,cmap=cmap)
			ax.axvline(x=10,ls='-',lw=1,c='w')
			if o == 0:
				ax.set_title(name,fontsize=10)
			if o == 2:
				xts = arange(-.1,.104,0.05)
				ax.axes.set_xticks([0,5,10,15,20])
				ax.axes.set_xticklabels(xts,fontsize=9)
				ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
				ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=7)
				if p != 0:
					ax.axes.set_yticklabels([])
			if p == 0:
				ax.set_ylabel('RMSD $\AA$',fontsize=9)
				yts = [int(round(i,roundery)) for i in midy][::int(float(len(midx))/ticknums)]
				ax.axes.set_yticks([[round(i,roundery) for i in midy].index(round(j,roundery)) for j in yts])
				ax.axes.set_yticklabels(yts,fontsize=7)
				if o != 2:
					ax.axes.set_xticklabels([])
			if o != 2 and p != 0:
				ax.axes.set_yticklabels([])
				ax.axes.set_xticklabels([])
	plt.savefig(pickles+'fig-dimple-sigma-vs-rmsd.png',dpi=500,bbox_inches='tight')
	plt.show()
	
#---compare the extents and residuals (2D histogram plots)
if do_errorlook_sigma:
	ticknums = 5
	rounderx = 2
	roundery = 0
	fig = plt.figure(figsize=(18,6))
	gs = gridspec.GridSpec(3,len(analyses),wspace=0.0,hspace=0.0)
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		order = (0,1,2)
		for o in order:
			ax = plt.subplot(gs[o,p])
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			target_zones = results_stack[p][order[o]].get(['type','target_zones'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
				for j in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][4])/10. < sigmamax]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][5])/10. < sigmamax]
			meansig = [1./2*(abs(params[i][4])/10.+abs(params[i][5])/10.) 
				for i in validhis if len(shape(params[i])) > 0]
			ax = plt.subplot(gs[o,p])
			H, xedges, yedges = histogram2d(meansig,resids,bins=21,range=((0,80),(0,20)),
				weights=[1./len(validhs) for i in validhs])
			midx = (xedges[1:]+xedges[:-1])/2
			midy = (yedges[1:]+yedges[:-1])/2
			extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
			cmap = mpl.cm.jet
			cmap.set_bad(cmap(0),1.)
			ax.imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
				norm=None,cmap=cmap)
			#ax.axvline(x=10,ls='-',lw=1,c='w')
			if o == 0:
				ax.set_title(name,fontsize=10)
			if o == 2:
				xts = arange(-.1,.104,0.05)
				xts = arange(0,81,20)
				ax.axes.set_xticks([0,5,10,15,20])
				ax.axes.set_xticklabels(arange(0,81,20),fontsize=9)
				ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
				ax.set_xlabel('$\mathsf{\sigma (nm)}$',fontsize=9)				
				if p != 0:
					ax.axes.set_yticklabels([])
			if p == 0:
				ax.set_ylabel('RMSD $\AA$',fontsize=9)
				yts = [int(round(i,roundery)) for i in midy][::int(float(len(midx))/ticknums)]
				ax.axes.set_yticks([[round(i,roundery) for i in midy].index(round(j,roundery)) for j in yts])
				ax.axes.set_yticklabels(yts,fontsize=7)
				if o != 2:
					ax.axes.set_xticklabels([])
			if o != 2 and p != 0:
				ax.axes.set_yticklabels([])
				ax.axes.set_xticklabels([])
	plt.savefig(pickles+'fig-dimple-sigma-vs-rmsd.png',dpi=500,bbox_inches='tight')
	plt.show()
	
#---compare the extents and curvatures (2D histogram plots)
if do_sigma_vs_hmax:
	ticknums = 5
	rounderx = 2
	roundery = 0
	fig = plt.figure(figsize=(18,6))
	gs = gridspec.GridSpec(3,len(analyses),wspace=0.0,hspace=0.0)
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		order = (0,1,2)
		for o in order:
			ax = plt.subplot(gs[o,p])
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			target_zones = results_stack[p][order[o]].get(['type','target_zones'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
				for j in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][4])/10. < sigmamax]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][5])/10. < sigmamax]
			meansig = [1./2*(abs(params[i][4])/10.+abs(params[i][5])/10.) 
				for i in validhis if len(shape(params[i])) > 0]
			ax = plt.subplot(gs[o,p])
			H, xedges, yedges = histogram2d(validhs,meansig,bins=21,range=(maxhrange,(0,80)),
				weights=[1./len(validhs) for i in validhs])
			midx = (xedges[1:]+xedges[:-1])/2
			midy = (yedges[1:]+yedges[:-1])/2
			extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
			cmap = mpl.cm.jet
			cmap.set_bad(cmap(0),1.)
			ax.imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
				norm=None,cmap=cmap)
			ax.axvline(x=10,ls='-',lw=1,c='w')
			if o == 0:
				ax.set_title(name,fontsize=10)
			if o == 2:
				xts = arange(-.1,.104,0.05)
				ax.axes.set_xticks([0,5,10,15,20])
				ax.axes.set_xticklabels(xts,fontsize=9)
				ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
				ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=7)
				if p != 0:
					ax.axes.set_yticklabels([])
			if p == 0:
				yts = [int(round(i,roundery)) for i in midy][::int(float(len(midx))/ticknums)]
				ax.axes.set_yticks([[round(i,roundery) for i in midy].index(round(j,roundery)) for j in yts])
				ax.axes.set_yticklabels(yts,fontsize=7)
				yts = arange(-.1,.104,0.05)
				yts = arange(0,81,20)
				ax.axes.set_yticks([0,5,10,15,20])
				ax.axes.set_yticklabels(arange(0,81,20),fontsize=9)
				ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
				ax.set_ylabel('$\mathsf{\sigma (nm)}$',fontsize=9)			
				if o != 2:
					ax.axes.set_xticklabels([])
			if o != 2 and p != 0:
				ax.axes.set_yticklabels([])
				ax.axes.set_xticklabels([])
	plt.savefig(pickles+'fig-dimple-sigma-vs-hmax.png',dpi=500,bbox_inches='tight')
	plt.show()
	
#---view the residual distributions
if do_resid_1d:
	ticknums = 5
	rounderx = 2
	roundery = 0
	minval = 0
	maxval = 15
	peakval = 0
	fig = plt.figure(figsize=(8,8))
	gs = gridspec.GridSpec(len(analyses),4,wspace=0.0,hspace=0.0)
	clrs2reord = [1,0,2]
	axes = []
	#---print report
	for p in range(len(analyses)):
		ccodes = analyses[p][1]
		name = analyses[p][2]
		print name
		fillcode = analyses[p][3]
		expected_direction = results_stack[p][0].notes[([i[0] 
			for i in results_stack[p][0].notes].index('expected_direction'))][1]
		order = ((0,1,2) if expected_direction == 1 else (1,0,2))
		order = (0,1,2)
		extraname = [' (-)','(+)','\n(unfiltered)']
		ax = plt.subplot(gs[p,0:4])
		for o in order:
			print 'expected_direction = '+str(expected_direction)
			print 'o = '+str(o)
			params = results_stack[p][order[o]].get(['type','params'])
			maxhs = results_stack[p][order[o]].get(['type','maxhs'])
			maxhxys = results_stack[p][order[o]].get(['type','maxhxys'])
			target_zones = results_stack[p][order[o]].get(['type','target_zones'])
			validhis = [i for i in range(len(maxhs)) 
				if (10*abs(maxhs[i]) > maxhfilter[0] and abs(10*maxhs[i]) < maxhfilter[1])]
			#---nanometer correction
			validhs = [10*maxhs[i] for i in validhis]
			resids = [sqrt(mean([abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
				for j in validhis]
			sigma_x = [abs(params[i][4])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][4])/10. < sigmamax]
			sigma_y = [abs(params[i][5])/10. for i in validhis if len(shape(params[i])) > 0 
				and abs(params[i][5])/10. < sigmamax]
			meansig = [1./2*(abs(params[i][4])/10.+abs(params[i][5])/10.) 
				for i in validhis if len(shape(params[i])) > 0]
			#---plots
			hist0,binedge0 = numpy.histogram(resids,bins=nbins,normed=False,
				weights=[1./len(validhs) for i in validhs],range=(minval,maxval))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			ax.plot(mid0,hist0,'o-',c=clrs2[clrs2reord[o]],alpha=1,lw=2,
				label=name+extraname[o])
			ax.legend(loc='upper right',prop={'size':10})
			ax.grid(True)
			axes.append(ax)
			if max(hist0) > peakval: peakval = max(hist0)
			ax.get_yaxis().set_major_locator(MaxNLocator(nbins=6,prune='both'))
			ax.axes.set_yticklabels([])
			if p != len(analyses)-1:
				ax.set_xticklabels([])
	for ax in axes:
		ax.set_ylim((0,1.1*peakval))
	axes[-1].set_xlabel('RMSD $\AA$',fontsize=14)
	plt.savefig(pickles+'fig-dimple-resids.png',dpi=500,bbox_inches='tight')
	plt.show()
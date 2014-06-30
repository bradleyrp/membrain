#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')
execfile('plotter.py')

import scipy.interpolate
import scipy.integrate
import subprocess

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import ndimage

'''
Script for reprocessing stress, framewise

Created 2013.12.03 by RPB

This program will open pickles with saved framewise first moment maps of the voxel-wise tension in bilayer 
	simulations and re-process them in order to study the distribution of spontaneous curvature. The previous
	script already incorporates bending rigidity, and it is important to remember that this study is only 
	qualitative at this point.

'''

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

exo70pip2_study = 1
enth_study = 0
#---Pickles containing kC0 plots from script-postproc-stress.py
if exo70pip2_study:
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
elif enth_study:        
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

#---plots
plot_maps = 1
plot_hist = 1
save_plots = 1
show_titles = 0

#---plot maps and histograms for +/- curvature domains (advanced version)
plot_plus_minus = 0

#---plot videos for +/- curvature domains (advanced version)
plot_plus_minus_video = 0
video_sharpfactor = 0.35
if plot_plus_minus_video:
	render_video_full = 1
	render_video_full_smooth = 1
	render_video_full_width = 40

#---file names
plot_advanced_maps_file = 'fig-'+outname+'-stress-framewise-maps-plusminus.png'
plot_advanced_hist_file = 'fig-'+outname+'-stress-framewise-histograms-plusminus.png'
plot_video_filebase = 'fig-'+outname+'-stress-framewise-maps-plusminus'
plot_maps_file = 'fig-'+outname+'-stress-framewise-maps.png'
plot_hist_file = 'fig-'+outname+'-stress-framewise-histograms.png'

#---Sign change
signchange = 1

#---LOAD
#-------------------------------------------------------------------------------------------------------------

raw_maps = []
for name in raw_maps_names:
	print 'Loading maps from '+name
	raw_maps.append(pickle.load(open(pickles+name,'r')))
msets = []
for name in pickle_structure_names:
	print 'Loading structures from '+name
	msets.append(pickle.load(open(pickles+name,'r')))
	
if save_plots: show_plots = 0

#---Plot maps / histograms for +/- curvature domains
#-------------------------------------------------------------------------------------------------------------

#---Note these are out of date - dropped after fixed errors in transpose, etc
if plot_plus_minus:	
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
	avg = [[],[],[]]
	clrs = brewer2mpl.get_map('Paired','qualitative', 8).mpl_colors
	#control_mean = mean([raw_maps[2][i][0] for i in range(len(raw_maps[2]))])+0.00048502604166666415
	#control_mean = mean([raw_maps[2][i][0] for i in range(len(raw_maps[2]))])
	control_mean = 0
	doms_sizes_poz = [[],[],[]]
	doms_sizes_neg = [[],[],[]]
	for k in range(3):
		print 'system = '+str(k)
		for i in range(len(raw_maps[k])):
			print 'frame = '+str(i)
			heatmap = raw_maps[k][i][0];
			tmp=ndimage.label(heatmap)
			#---use this and it looks like control is smoother on histogram
			#im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*40))
			#---use this and it looks opposite
			#im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*50))
			#---another filtering attempt
			im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*20),mode='wrap')
			maskp = im > im.mean()
			maskp = im > control_mean
			doms_sizes_poz[k].append(sum(im > im.mean())/4096.)
			doms_sizes_neg[k].append(sum(im < im.mean())/4096.)
			avg[k].append(maskp)
	fig = plt.figure()
	numgridpts = 64
	avgmeans = [mean(avg[i],axis=0)-0.5 for i in range(3)]
	#avgmeans = [std(avg[i],axis=0) for i in range(3)]
	#avgmeans = [avg[i][0]-0.5 for i in range(3)]
	extremum = max([max([max(abs(i)) for i in avgmeans[j]]) for j in range(3)])

	#---Plot average +/- map
	ax0 = plt.subplot2grid((1,4), (0,0))
	ax0.set_title(mapslabels[0])
	'''
	protein_centers = result_stack[0][1]
	nprots = nprots_list[0]
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax0.add_patch(circ)
	'''
	protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in np.mean(msets[0].protein,axis=0)])
	hull = scipy.spatial.ConvexHull(protpts)
	ax0.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'r-',lw=1)
	ax0.plot(protpts[hull.vertices,0][-1],protpts[hull.vertices,1][0],'r-',lw=1)
	#---end edit
	ax0.imshow(array(avgmeans[0]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax1 = plt.subplot2grid((1,4), (0,1))
	ax1.set_title(mapslabels[1])
	'''
	protein_centers = result_stack[1][1]
	nprots = nprots_list[1]
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax1.add_patch(circ)
	'''
	protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] 
		for i in np.mean(msets[1].protein,axis=0)])
	hull = scipy.spatial.ConvexHull(protpts)
	ax1.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'r-',lw=1)
	ax1.plot(protpts[hull.vertices,0][-1],protpts[hull.vertices,1][0],'r-',lw=1)
	#---end edit
	ax1.imshow(array(avgmeans[1]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax2 = plt.subplot2grid((1,4), (0,2))
	ax2.set_title(mapslabels[2])
	img = ax2.imshow(array(avgmeans[2]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,
		vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	cax = inset_axes(ax2,
         width="5%",
         height="100%",
         bbox_transform=ax2.transAxes,
         bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
         loc= 1)
	fig.colorbar(img,cax=cax)
	cax.tick_params(labelsize=10) 
	cax.set_ylabel('occupancy',fontsize=10)
	plt.tight_layout() 
	if show_plots:
		plt.show()
	if save_plots:
		plt.savefig(pickles+plot_advanced_maps_file,dpi=500,bbox_inches='tight')
	plt.cla()

	#---Plot +/- histograms
	fig = plt.figure()
	plt.grid(True)
	nbins = 10
	minsize = 0
	plt.xlabel(r'fraction positive spontaneous curvature',labelpad = 10,fontsize=20)
	plt.ylabel('frequency', labelpad = 10,fontsize=20)
	hist0 = numpy.histogram(avgmeans[0],bins=nbins)
	hist1 = numpy.histogram(avgmeans[1],bins=nbins)
	hist2 = numpy.histogram(avgmeans[2],bins=nbins)
	plt.plot(hist0[1][1:],hist0[0],color=clrs[1],alpha=1.,lw=2,
		label=mapslabels[0])
	plt.plot(hist1[1][1:],hist1[0],color=clrs[3],alpha=1.,lw=2,
		label=mapslabels[1])
	plt.plot(hist2[1][1:],hist2[0],color=clrs[5],alpha=1.,lw=2,
		label=mapslabels[2])
	plt.legend()
	plt.tight_layout() 
	if show_plots:
		plt.show()
	if save_plots:
		plt.savefig(pickles+plot_advanced_hist_file,dpi=500,bbox_inches='tight')
	plt.cla()
		
#---Plot videos for +/- curvature domains
#-------------------------------------------------------------------------------------------------------------

if plot_plus_minus_video:	
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
	avg = [[],[],[]]
	clrs = brewer2mpl.get_map('Paired', 'qualitative', 8).mpl_colors
	control_mean = 0
	doms_sizes_poz = [[],[],[]]
	doms_sizes_neg = [[],[],[]]
	for k in range(3):
		print 'system = '+str(k)
		for i in range(len(raw_maps[k])):
			print 'frame = '+str(i)
			heatmap = raw_maps[k][i][0];
			tmp=ndimage.label(heatmap)
			#---use this and it looks like control is smoother on histogram
			#im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*40))
			#---use this and it looks opposite
			#im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*50))
			#---another filtering attempt
			im = ndimage.gaussian_filter(heatmap, sigma=256/(4.*20),mode='wrap')
			maskp = im > im.mean()
			maskp = im > control_mean
			doms_sizes_poz[k].append(sum(im > im.mean())/4096.)
			doms_sizes_neg[k].append(sum(im < im.mean())/4096.)
			avg[k].append(maskp)
	print shape(avg[0])
	print shape(avg[1])
	print shape(avg[2])
	print 'bkah'
	avgmeans_stack = list(array([[avg[i][j]-0.5 for j in range(len(avg[i]))] for i in range(3)]).T)
	print 'precomputing maps ...'
	if render_video_full:
		precomp_maps = list(array([[raw_maps[k][fr][0] for fr in range(len(avg[k]))] for k in range(3)]).T)
		extremum = 0.5
	else:
		precomp_maps = avgmeans_stack[fr]
		extremum = 0.5
	if render_video_full_smooth:
		precomp_maps = array([[ndimage.gaussian_filter(precomp_maps[k][fr],
			sigma=256/(4.*render_video_full_width)) for fr in range(len(avg[k]))] for k in range(3)]).T
		extremum = max([np.max([np.max(i) for i in precomp_maps]),np.abs(np.max([np.max(i) 
			for i in precomp_maps]))])
		#---take the half extremum to make it look sharper
		extremum = video_sharpfactor*extremum
		print 'extremum = '+str(extremum)
	for fr in range(min([len(avg[i]) for i in range(3)])):
		print 'rendering frame '+str(fr)
		#avgmeans = precomp_maps[fr]
		avgmeans = [precomp_maps[i][fr] for i in range(3)]
		fig = plt.figure()
		numgridpts = 64
		ax0 = plt.subplot2grid((1,3), (0,0))
		ax0.set_title(mapslabels[0])
		'''
		nprots = nprots_list[0]
		if nprots > 0:
			protlen = int(shape(msets[0].protein[0])[0]/nprots)
			protlenshow = protlen
			protein_centers = mean([msets[0].protein[fr][i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)],axis=1)
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax0.add_patch(circ)
		'''
		protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in msets[0].protein[fr]])
		hull = scipy.spatial.ConvexHull(protpts)
		ax0.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'r-',lw=1)
		ax0.plot(protpts[hull.vertices,0][-1],protpts[hull.vertices,1][0],'r-',lw=1)
		#---end edit
		ax0.imshow(array(avgmeans[0]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		ax1 = plt.subplot2grid((1,3), (0,1))
		ax1.set_title(mapslabels[1])
		protein_centers = result_stack[1][1]
		'''
		nprots = nprots_list[1]
		if nprots > 0:
			protlen = int(shape(msets[1].protein[0])[0]/nprots)
			protlenshow = protlen
			protein_centers = mean([msets[1].protein[fr][i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)],axis=1)
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax1.add_patch(circ)
		'''
		protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in msets[1].protein[fr]])
		hull = scipy.spatial.ConvexHull(protpts)
		ax1.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'r-',lw=1)
		ax1.plot(protpts[hull.vertices,0][-1],protpts[hull.vertices,1][0],'r-',lw=1)
		#---end edit
		ax1.imshow(array(avgmeans[1]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		ax2 = plt.subplot2grid((1,3), (0,2))
		ax2.set_title(mapslabels[2])
		img = ax2.imshow(array(avgmeans[2]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		plt.tight_layout() 
		plt.savefig(pickles+plot_video_filebase+'.fr.'+str('%04d'%fr)+'.png',dpi=500,bbox_inches='tight')
		plt.cla()
		plt.close()
	print 'done rendering. run this:'
	print r"ffmpeg -i stress-framewise-plusminus-maps.fr.%04d.png -vcodec mpeg1video -qscale 0  ",
	print r"-filter:v 'setpts=2.0*PTS' ./untitled.mpg"
	subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png','-vcodec','mpeg1video',
		'-qscale','0','-filter:v','setpts=2.0*PTS',pickles+'/'+plot_video_filebase+'.mpeg'])
	for fr in range(min([len(avg[i]) for i in range(3)])): 
		subprocess.call(['rm',pickles+'/'+plot_video_filebase+'.fr.'+str('%04d'%fr)+'.png'])


#---Plot spontaneous curvature (C0) histograms, averaged across all voxels and all frames (together)
#-------------------------------------------------------------------------------------------------------------

if plot_hist:
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
	plt.plot(mid0,hist0,'bo-',alpha=1.,lw=2,label=mapslabels[0])
	plt.plot(mid1,hist1,'co-',alpha=1.,lw=2,label=mapslabels[1])
	plt.plot(mid2,hist2,'ko-',alpha=1.,lw=2,label=mapslabels[2])
	plt.xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',labelpad = 10,fontsize=20)
	plt.ylabel('frequency', labelpad = 10,fontsize=20)
	ax.set_xlim((-0.05,0.05))
	if show_titles:
		plt.title('spontaneous curvature')
	fig.tight_layout()
	plt.legend()
	plt.tight_layout() 
	if save_plots:
		plt.savefig(pickles+plot_hist_file,dpi=500,bbox_inches='tight')
	if show_plots:
		plt.show()
	plt.cla()
	plt.close()

#---Plot spontaneous curvature (C0) maps, averaged over all frames
#-------------------------------------------------------------------------------------------------------------

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
	extremum = max([max([max(i) for i in result_stack[j][0]]) for j in range(3)])
	extremum = 0.06
	ax0 = plt.subplot2grid((1,4),(0,0))
	ax0.set_title(mapslabels[0])
	ax0.set_xticklabels([])
	ax0.set_yticklabels([])
	ax0.set_adjustable('box-forced')
	dat = result_stack[0][0]
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
			ax0.plot(shifthullx,shifthully,'k-',lw=0.6)
	ax0.imshow(signchange*array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax1 = plt.subplot2grid((1,4),(0,1))
	ax1.set_title(mapslabels[1])
	ax1.set_xticklabels([])
	ax1.set_yticklabels([])
	ax1.set_adjustable('box-forced')
	dat = result_stack[1][0]
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
			ax1.plot(shifthullx,shifthully,'k-',lw=0.6)
	ax1.imshow(signchange*array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax2 = plt.subplot2grid((1,4), (0,2),colspan=1)
	ax2.set_title(mapslabels[2])
	ax2.set_xticklabels([])
	ax2.set_yticklabels([])
	ax2.set_adjustable('box-forced')
	dat = result_stack[2][0]
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
	if save_plots:
		plt.savefig(pickles+plot_maps_file,dpi=500,bbox_inches='tight')
	if show_plots:
		plt.show()
	plt.cla()
	plt.close()


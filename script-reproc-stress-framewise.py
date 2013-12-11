#!/usr/bin/python -i

from membrainrunner import *

location = 'light'
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

Needs: video
'''

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Pickles containing kC0 plots from script-postproc-stress.py
raw_maps_names = [
	'pkl.stressdecomp.membrane-v614.part0002.pkl',
	'pkl.stressdecomp.membrane-v612.part0003.pkl',
	'pkl.stressdecomp.membrane-v550.part0008.pkl']
pickle_structure_names = [
	'pkl.structures.membrane-v614-stress.md.part0002.pkl',
	'pkl.structures.membrane-v612-stress.md.part0003.pkl',
	'pkl.structures.membrane-v550-stress.md.part0008.shifted.pkl']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

raw_maps = []
for name in raw_maps_names:
	print 'Loading maps from '+name
	raw_maps.append(pickle.load(open(pickles+name,'r')))
msets = []
for name in pickle_structure_names:
	print 'Loading structures from '+name
	msets.append(pickle.load(open(pickles+name,'r')))

#---plots
plot_maps = 0
plot_hist = 0
#---plot the domain size for connected positive and negative curvature domains (version 1)
plot_connected_domains_basic = 0
#---plot the domain size for connected positive and negative curvature domains (version 2)
plot_plus_minus_basic = 0
if plot_plus_minus_basic:
	plot_plus_minus_stacked = 1
	plot_plus_minus_both = 1
#---plot maps, histograms, and videos for +/- curvature domains (advanced version)
plot_plus_minus = 0
if plot_plus_minus:
	plot_advanced_hist0 = 0
	plot_advanced_hist = 1
plot_plus_minus_video = 1
if plot_plus_minus_video:
	render_video_full = 1
	render_video_full_smooth = 1
	render_video_full_width = 5

#---plotting parameters
save_plots = 1
if save_plots:
	show_plots = 0
show_titles = 0

#---file names
plot_maps_file = 'stress-framewise-comparison.png'
plot_hist_file = 'stress-framewise-histograms.png'
plot_connected_domains_basic_file = 'stress-framewise-domains.png'
plot_plus_minus_stacked_file = 'stress-framewise-histograms-plus-minus.png'
plot_plus_minus_both_file = 'stress-framewise-histograms-plus-minus.png'
plot_advanced_maps_file = 'stress-framewise-maps.png'
plot_advanced_hist_file = 'stress-framewise-histograms.png'
plot_video_filebase = 'stress-framewise-plusminus-maps'

#---plot maps, histograms, and videos for +/- curvature domains (advanced version)
if plot_plus_minus:	
	nprots_list = [4,2,0]
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
	control_mean = mean([raw_maps[2][i][0] for i in range(len(raw_maps[2]))])+0.00048502604166666415
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
	ax0 = plt.subplot2grid((1,4), (0,0))
	ax0.set_title(r'$\textbf{{ENTH}\ensuremath{\times}4}$')
	protein_centers = result_stack[0][1]
	nprots = 4
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax0.add_patch(circ)
	ax0.imshow(array(avgmeans[0]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax1 = plt.subplot2grid((1,4), (0,1))
	ax1.set_title(r'$\textbf{{ENTH}\ensuremath{\times}1}$')
	protein_centers = result_stack[1][1]
	nprots = 2
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax1.add_patch(circ)
	ax1.imshow(array(avgmeans[1]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax2 = plt.subplot2grid((1,4), (0,2))
	ax2.set_title(r'$\textbf{{control}}$')
	img = ax2.imshow(array(avgmeans[2]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
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
	if plot_advanced_hist0:
		nbins = 15
		minsize = 0
		hist0 = numpy.histogram([i for i in concatenate((doms_sizes_neg[0],doms_sizes_poz[0])) 
			if i > minsize],bins=nbins)
		hist1 = numpy.histogram([i for i in concatenate((doms_sizes_neg[1],doms_sizes_poz[1])) 
			if i > minsize],bins=nbins)
		hist2 = numpy.histogram([i for i in concatenate((doms_sizes_neg[2],doms_sizes_poz[2])) 
			if i > minsize],bins=nbins)
		plt.plot(hist0[1][1:],hist0[0],color=clrs[1],alpha=1.,lw=2,
			label=r'$\textbf{+ {ENTH}\ensuremath{\times}4}$')
		plt.plot(hist1[1][1:],hist1[0],color=clrs[3],alpha=1.,lw=2,
			label=r'$\textbf{+ {ENTH}\ensuremath{\times}1}$')
		plt.plot(hist2[1][1:],hist2[0],color=clrs[5],alpha=1.,lw=2,
			label=r'$\textbf{+ {control}}$')
		plt.legend()
		plt.tight_layout() 
		if show_plots:
			plt.show()
		#if save_plots:
		#	plt.savefig(pickles+plot_maps_file_advanced,dpi=500,bbox_inches='tight')
		plt.cla()
	if plot_advanced_hist:
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
			label=r'$\textbf{+ {ENTH}\ensuremath{\times}4}$')
		plt.plot(hist1[1][1:],hist1[0],color=clrs[3],alpha=1.,lw=2,
			label=r'$\textbf{+ {ENTH}\ensuremath{\times}1}$')
		plt.plot(hist2[1][1:],hist2[0],color=clrs[5],alpha=1.,lw=2,
			label=r'$\textbf{+ {control}}$')
		plt.legend()
		plt.tight_layout() 
		if show_plots:
			plt.show()
		if save_plots:
			plt.savefig(pickles+plot_advanced_hist_file,dpi=500,bbox_inches='tight')
		plt.cla()
		
#---plot maps, histograms, and videos for +/- curvature domains (advanced version)
if plot_plus_minus_video:	
	nprots_list = [4,2,0]
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
	control_mean = mean([raw_maps[2][i][0] for i in range(len(raw_maps[2]))])+0.00048502604166666415
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
	avgmeans_stack = [[avg[i][j]-0.5 for i in range(3)] for j in range(len(avg[0]))]
	print 'precomputing maps ...'
	if render_video_full:
		precomp_maps = [[raw_maps[k][fr][0] for k in range(3)] for fr in range(len(avg[0]))]
		extremum = 0.5
	else:
		precomp_maps = avgmeans_stack[fr]
		extremum = 0.5
	if render_video_full_smooth:
		precomp_maps = [[ndimage.gaussian_filter(precomp_maps[fr][i],
			sigma=256/(4.*render_video_full_width)) for i in range(3)] for fr in range(len(avg[0]))]
		extremum = np.max(np.abs(precomp_maps))
		print 'extremum = '+str(extremum)
	for fr in range(len(avg[0])):
		print 'rendering frame '+str(fr)
		avgmeans = precomp_maps[fr]
		fig = plt.figure()
		numgridpts = 64
		ax0 = plt.subplot2grid((1,3), (0,0))
		ax0.set_title(r'$\textbf{{ENTH}\ensuremath{\times}4}$')
		nprots = 4
		if nprots > 0:
			protlen = int(shape(msets[0].protein[0])[0]/nprots)
			protein_centers = mean([msets[0].protein[fr][i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)],axis=1)
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax0.add_patch(circ)
		ax0.imshow(array(avgmeans[0]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		ax1 = plt.subplot2grid((1,3), (0,1))
		ax1.set_title(r'$\textbf{{ENTH}\ensuremath{\times}1}$')
		protein_centers = result_stack[1][1]
		nprots = 1
		if nprots > 0:
			protlen = int(shape(msets[1].protein[0])[0]/nprots)
			protein_centers = mean([msets[1].protein[fr][i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)],axis=1)
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax1.add_patch(circ)
		ax1.imshow(array(avgmeans[1]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		ax2 = plt.subplot2grid((1,3), (0,2))
		ax2.set_title(r'$\textbf{{control}}$')
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
	subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png','-vcodec','mpeg1video','-qscale','0','-filter:v','setpts=2.0*PTS',pickles+'/'+plot_video_filebase+'.mpeg'])
	for fr in range(len(avg[0])): subprocess.call(['rm',pickles+'/'+plot_video_filebase+'.fr.'+str('%04d'%fr)+'.png'])

#---Plot the domain size for connected positive and negative curvature domains (version 1)
if plot_connected_domains_basic:
	doms_sizes = [[],[],[]]
	for k in range(3):
		print 'system = '+str(k)
		for i in range(len(raw_maps[k])):
			print 'frame = '+str(i)
			heatmap = raw_maps[k][i][0];
			disc = array([[(1 if heatmap[i][j]>0 else 0) 
				for i in range(shape(heatmap)[0])] 
				for j in range(shape(heatmap)[1])])
			doms,num = ndimage.measurements.label(1-disc)
			doms_sizes[k].extend([sum(doms==i) for i in range(num)])
	nbins = 12
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.rc('font', family='sans-serif')
	ax.grid(True)
	hist0 = numpy.histogram([i for i in doms_sizes[0] if i > 1500],bins=nbins)
	hist1 = numpy.histogram([i for i in doms_sizes[1] if i > 1500],bins=nbins)
	hist2 = numpy.histogram([i for i in doms_sizes[2] if i > 1500],bins=nbins)
	plt.plot(hist0[1][1:],hist0[0],color='b',alpha=1.,lw=2,label=r'$\textbf{{ENTH}\ensuremath{\times}4}$')
	plt.plot(hist1[1][1:],hist0[0],color='c',alpha=1.,lw=2,label=r'$\textbf{{ENTH}\ensuremath{\times}1}$')
	plt.plot(hist2[1][1:],hist0[0],color='k',alpha=1.,lw=2,label=r'$\textbf{{control}}$')
	ax.set_ylim((0,50))
	plt.xlabel(r'connected domain size $\mathsf{(nm^{2})}$',labelpad = 10,fontsize=20)
	plt.ylabel('frequency', labelpad = 10,fontsize=20)
	if show_titles:
		plt.title(r'connected $+/-$ curvature domains')
	plt.legend()
	plt.tight_layout() 
	if save_plots:
		plt.savefig(pickles+plot_connected_domains_basic_file,dpi=500,bbox_inches='tight')
	if show_plots:
		plt.show()
	plt.cla()
	plt.close()
		
#---Basic method for computing positive vs negative spontaneous curvature (C0) domain sizes		
if plot_plus_minus_basic:		
	thresh = 0.02
	doms_sizes_neg = [[],[],[]]
	for k in range(3):
		print 'system = '+str(k)
		for i in range(len(raw_maps[k])):
			print 'frame = '+str(i)
			heatmap = raw_maps[k][i][0];
			disc0 = array([[(1 if abs(heatmap[i][j]) > thresh else 0) 
				for i in range(shape(heatmap)[0])] 
				for j in range(shape(heatmap)[1])])
			disc1 = array([[(1 if heatmap[i][j] > 0 else -1) 
				for i in range(shape(heatmap)[0])] 
				for j in range(shape(heatmap)[1])])
			disc3 = disc0 + disc1
			doms,num = ndimage.measurements.label(disc3==0)
			doms_sizes_neg[k].extend([sum(doms==i) for i in range(num)])
	doms_sizes_poz = [[],[],[]]
	for k in range(3):
		print 'system = '+str(k)
		for i in range(len(raw_maps[k])):
			print 'frame = '+str(i)
			heatmap = raw_maps[k][i][0];
			disc0 = array([[(1 if abs(heatmap[i][j]) > thresh else 0) 
				for i in range(shape(heatmap)[0])] 
				for j in range(shape(heatmap)[1])])
			disc1 = array([[(1 if heatmap[i][j] > 0 else -1) 
				for i in range(shape(heatmap)[0])] 
				for j in range(shape(heatmap)[1])])
			disc3 = disc0 + disc1
			doms,num = ndimage.measurements.label(disc3==2)
			doms_sizes_poz[k].extend([sum(doms==i) for i in range(num)])
	#---Plots stacked histograms with positive and negative domain sizes for each system
	clrs = brewer2mpl.get_map('Paired', 'qualitative', 8).mpl_colors
	nbins = 10
	minsize = 1500
	hist0a = numpy.histogram([i for i in doms_sizes_neg[0] if i > minsize],bins=nbins)
	hist1a = numpy.histogram([i for i in doms_sizes_neg[1] if i > minsize],bins=nbins)
	hist2a = numpy.histogram([i for i in doms_sizes_neg[2] if i > minsize],bins=nbins)
	hist0b = numpy.histogram([i for i in doms_sizes_poz[0] if i > minsize],bins=nbins)
	hist1b = numpy.histogram([i for i in doms_sizes_poz[1] if i > minsize],bins=nbins)
	hist2b = numpy.histogram([i for i in doms_sizes_poz[2] if i > minsize],bins=nbins)
	fig = plt.figure()
	plt.title(r'continuous peak domain area')
	ax0 = plt.subplot2grid((3,1),(0,0))
	ax0.plot(hist0a[1][1:],hist0a[0],color=clrs[1],alpha=1.,lw=2,
		label=r'$\textbf{+ {ENTH}\ensuremath{\times}4}$')
	ax0.plot(hist0b[1][1:],hist0b[0],color=clrs[0],alpha=1.,lw=2,
		label=r'$\textbf{- {ENTH}\ensuremath{\times}4}$')
	ax0.legend(loc=2)
	ax0.set_ylim((0,50))
	#ax0.set_xlim((3500,3900))
	ax1 = plt.subplot2grid((3,1),(1,0))
	ax1.plot(hist1a[1][1:],hist1a[0],color=clrs[3],alpha=1.,lw=2,
		label=r'$\textbf{+ {ENTH}\ensuremath{\times}1}$')
	ax1.plot(hist1b[1][1:],hist1b[0],color=clrs[2],alpha=1.,lw=2,
		label=r'$\textbf{- {ENTH}\ensuremath{\times}1}$')
	ax1.legend(loc=2)
	ax1.set_ylim((0,50))
	#ax1.set_xlim((3500,3900))
	ax2 = plt.subplot2grid((3,1),(2,0))
	ax2.plot(hist2a[1][1:],hist2a[0],color=clrs[5],alpha=1.,lw=2,label=r'$\textbf{+ {control}}$')
	ax2.plot(hist2b[1][1:],hist2b[0],color=clrs[4],alpha=1.,lw=2,label=r'$\textbf{- {control}}$')
	ax2.legend(loc=2)
	ax2.set_ylim((0,50))
	#ax2.set_xlim((3500,3900))
	plt.xlabel(r'domain area $(\textbf{{nm}^{2}})$',labelpad = 10,fontsize=16)
	plt.tight_layout() 
	if show_plots:
		plt.show()
	if save_plots:
		plt.savefig(pickles+plot_plus_minus_stacked_file,dpi=500,bbox_inches='tight')
	plt.cla()
	plt.close()
	#---Plots histograms for the total domain size for each system
	nbins = 10
	minsize = 1500
	hist0 = numpy.histogram([i for i in concatenate((doms_sizes_neg[0],doms_sizes_poz[0])) 
		if i > minsize],bins=nbins)
	hist1 = numpy.histogram([i for i in concatenate((doms_sizes_neg[1],doms_sizes_poz[1])) 
		if i > minsize],bins=nbins)
	hist2 = numpy.histogram([i for i in concatenate((doms_sizes_neg[2],doms_sizes_poz[2])) 
		if i > minsize],bins=nbins)
	plt.plot(hist0[1][1:],hist0[0],color=clrs[1],alpha=0.8,lw=2,
		label=r'$\textbf{+ {ENTH}\ensuremath{\times}4}$')
	plt.plot(hist1[1][1:],hist1[0],color=clrs[3],alpha=0.8,lw=2,
		label=r'$\textbf{+ {ENTH}\ensuremath{\times}1}$')
	plt.plot(hist2[1][1:],hist2[0],color=clrs[5],alpha=0.8,lw=2,
		label=r'$\textbf{+ {control}}$')
	plt.xlabel(r'domain area $(\textbf{{nm}^{2}})$',labelpad = 10,fontsize=16)
	plt.ylabel('frequency', labelpad = 10,fontsize=20)
	plt.legend()
	plt.grid(True)
	plt.tight_layout() 
	if show_plots:
		plt.show()
	if save_plots:
		plt.savefig(pickles+plot_plus_minus_both_file,dpi=500,bbox_inches='tight')
	plt.cla()
	plt.close()

#---Plot spontaneous curvature (C0) distribution, averaged across all voxels and all frames (together)
if plot_hist:
	pdist0 = [i for j in mean([i[0] for i in raw_maps[0]],axis=0) for i in j]
	pdist1 = [i for j in mean([i[0] for i in raw_maps[1]],axis=0) for i in j]
	pdist2 = [i for j in mean([i[0] for i in raw_maps[2]],axis=0) for i in j]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.rc('font', family='sans-serif')
	ax.grid(True)
	nbins = 20
	hist0 = numpy.histogram(pdist0,bins=nbins)
	hist1 = numpy.histogram(pdist1,bins=nbins)
	hist2 = numpy.histogram(pdist2,bins=nbins)
	plt.plot(hist0[1][1:],hist0[0],color='b',alpha=1.,lw=2,label=r'$\textbf{{ENTH}\ensuremath{\times}4}$')
	plt.plot(hist1[1][1:],hist1[0],color='c',alpha=1.,lw=2,label=r'$\textbf{{ENTH}\ensuremath{\times}1}$')
	plt.plot(hist2[1][1:],hist2[0],color='k',alpha=1.,lw=2,label=r'$\textbf{{control}}$')
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
if plot_maps:
	nprots_list = [4,2,0]
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
	ax0 = plt.subplot2grid((1,4),(0,0))
	ax0.set_title(r'$\textbf{{ENTH}\ensuremath{\times}4}$')
	ax0.set_xticklabels([])
	ax0.set_yticklabels([])
	ax0.set_adjustable('box-forced')
	dat = result_stack[0][0]
	protein_centers = result_stack[0][1]
	nprots = 4
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax0.add_patch(circ)
	ax0.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax1 = plt.subplot2grid((1,4),(0,1))
	ax1.set_title(r'$\textbf{{ENTH}\ensuremath{\times}1}$')
	ax1.set_xticklabels([])
	ax1.set_yticklabels([])
	ax1.set_adjustable('box-forced')
	dat = result_stack[1][0]
	protein_centers = result_stack[1][1]
	nprots = 2
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax1.add_patch(circ)
	ax1.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	ax2 = plt.subplot2grid((1,4), (0,2),colspan=1)
	ax2.set_title(r'$\textbf{{control}}$')
	ax2.set_xticklabels([])
	ax2.set_yticklabels([])
	ax2.set_adjustable('box-forced')
	dat = result_stack[2][0]
	img = ax2.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
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


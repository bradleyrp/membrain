#!/usr/bin/python -i

from membrainrunner import *

import numpy as N
import pylab
from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq
import os

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---settings
skip = 1
framecount = None
location = ''
execfile('locations.py')

#---parameters
height_direction = 1
cutoff_distance = 5.
cutoff_distance = 15.
make_figs = False

#---analysis plan
analysis_plan = slice(None,None)
analysis_descriptors = [
	('pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',slice(None),None,''),
	('pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',slice(None),None,''),
	('pkl.structures.membrane-v612-stress.md.part0003.pkl',slice(None),None,''),
	('pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',slice(None),None,''),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
	'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl','.prot-v614'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
	'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl','.prot-v700')]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def blockplot(blocks,zvals=None,colorcodes=None,filename=None,show=True):
	'''
	BLOCKPLOT
	Takes a 2D array and plots it on a grid.
	Options:
	zvals : a dictionary which maps values of blocks to an integer range
	colorcodes : map from an integer range to RGB colors
	blocks: a 2D array of values
	'''
	if colorcodes != None and zvals == None:
		zvals = range(len(colors))
	dims = np.shape(blocks)
	fig = plt.figure()
	ax = fig.add_subplot(111,aspect='equal')
	ax.set_xlim((0,dims[0]))
	ax.set_ylim((0,dims[1]))
	colors = 100*np.random.random(dims[0]*dims[1])
	patches = []
	plotcolors = []
	for x in range(dims[0]):
		for y in range(dims[1]):
			if colorcodes == None:
				rect = mpl.patches.Rectangle((x,y),1,1,linewidth=2.0,linestyle='none')
			else:
				rect = mpl.patches.Rectangle((x,y),1,1,linewidth=0.0,facecolor=colorcodes[zvals[blocks[x,y]]])
			patches.append(rect)
			plotcolors.append(blocks[x,y])
	if colorcodes == None:
		p = mpl.collections.PatchCollection(patches,
			cmap=brewer2mpl.get_map('Paired','qualitative',4,'reverse').mpl_colormap,
			alpha=1.0,edgecolors='none')
		p.set_array(array(plotcolors))
	else:
		p = mpl.collections.PatchCollection(patches,alpha=1.0,match_original=True)
	ax.add_collection(p)
	if filename != None:
		plt.savefig(filename)
	if show:
		plt.show()
	else:
		plt.close()
	
def lateral_discretize(fr,result='area'):
	'''
	LATERAL_DISCRETIZE
	Discretize a boolean of bilayer midplane sign given access to a membrane object and a proteins object.
	Return protein present/absent region.
	Options:
	result : default return just the area calculation, grid returns the result grids, all returns both
	'''
	print 'Analyzing discretized areas for frame: '+str(fr)
	#---transpose correction
	surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] 
		for i in range(mset.griddims[0]-1)]).T
	protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
		dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
	if result == 'all':
		return [[float(sum(surf_discrete==1)),float(sum(surf_discrete==-1)),
			float(sum(surf_discrete+buf==2)),float(sum(surf_discrete+buf==0)),float(sum(buf==1))],
			surf_discrete,buf]
	elif result == 'grid':
		return [surf_discrete,buf]
	else:
		return [float(sum(surf_discrete==1)),float(sum(surf_discrete==-1)),
			float(sum(surf_discrete+buf==2)),float(sum(surf_discrete+buf==0)),float(sum(buf==1))]
			
def view_example(fr=0):
	'''Diagnostic plot of frame filter figure.'''
	result = lateral_discretize(fr,result='grid')
	blockplot(result[0]+result[1],zvals=zvals,colorcodes=colorcodes)

def view_example_mean(filename=None):
	'''Diagnostic plot of frame filter figure, mean configuration.'''
	mset.calculate_average_surface()
	surf_discrete = array([[(1 if mset.surf_mean[i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] 
		for i in range(mset.griddims[0]-1)])
	protpos = array(where(array(mset.rezipgrid(mean(proteins,axis=0)))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
		dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
	blockplot(surf_discrete+buf,zvals=zvals,colorcodes=colorcodes,filename=filename,show=False)
	
def view_figures(area_counts):
	'''Generate figures and print results.'''
	mpl.rcParams.update({'font.size': 30})
	area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((2,1),(0,0))
	t = range(len(area_counts))
	posarea = array([area_per_tile*area_counts[i,2] for i in range(len(area_counts))])
	negarea = array([area_per_tile*area_counts[i,3] for i in range(len(area_counts))])
	ax.plot(posarea,'r-',label='$z>0$',lw=3)
	ax.plot(negarea,'b-',label='$z<0$',lw=3)
	ax.fill_between(t, negarea,posarea, facecolor=colorcodes[zvals[-1]], alpha=0.5,where=negarea>posarea)
	ax.fill_between(t, posarea,negarea, facecolor=colorcodes[zvals[1]], alpha=0.5,where=negarea<posarea)
	ax.set_xlim((0,len(area_counts)))
	#ax.spines['top'].set_visible(False)
	#ax.spines['right'].set_visible(False)
	#ax.spines['bottom'].set_position(('outward', 20))
	#ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 	
	plt.xlabel('Frame', labelpad = 10)
	ylabel1 =plt.ylabel('Area ($nm^2$)', labelpad = 10)
	#plt.title('Tile area (within $'+str(cutoff_distance)+'$ nm of protein) by midplane height')
	legend1 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax.grid()
	ax2 = plt.subplot2grid((2,1),(1,0))
	histplot2 = ax2.hist(posarea,color=colorcodes[zvals[1]],alpha=0.8,label='$z>0$')
	histplot1 = ax2.hist(negarea,color=colorcodes[zvals[-1]],alpha=0.8,label='$z<0$')
	#ax2.spines['top'].set_visible(False)
	#ax2.spines['right'].set_visible(False)
	#ax2.spines['bottom'].set_position(('outward', 20))
	#ax2.spines['left'].set_position(('outward', 30))
	for tick in ax2.xaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	for tick in ax2.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 	
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ylabel2 = plt.ylabel('Frames', labelpad = 10)
	plt.xlabel('Area ($nm^2$)', labelpad = 10)
	legend2 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax2.grid()
	plt.tight_layout()
	#plt.savefig(pickles+'result.fig.tilefilter.areas.'+systemprefix+'.png',dpi=300,
	#	bbox_extra_artists=[ylabel1,ylabel2,legend1,legend2],bbox_inches='tight')
	plt.savefig(pickles+'fig-'+sysname+'-tilefilter.areas.png',dpi=300,
		bbox_extra_artists=[ylabel1,ylabel2,legend1,legend2],bbox_inches='tight')
	plt.close()
	#fp = open(pickles+'result.txt.tilefilter.'+systemprefix+'.txt','w')
	fp = open(pickles+'dat-'+sysname+'-tilefilter.areas.txt','w')
	fp.write('mean total area: '+str(product(vecs[0:2])/100.)+'nm2\n')
	fp.write('mean positive area: '+str(mean(area_per_tile*area_counts[:,0]))+' nm2\n')
	fp.write('mean negative area: '+str(mean(area_per_tile*area_counts[:,1]))+' nm2\n')
	fp.write('mean protein/positive area: '+str(area_per_tile*mean(area_counts[:,2]))+' nm2\n')
	fp.write('mean protein/negative area: '+str(area_per_tile*mean(area_counts[:,3]))+' nm2\n')
	fp.close()
	
def view_figures_sweep(area_counts_sweep):	
	mpl.rcParams.update({'font.size': 30})
	area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((1,1),(0,0))
	if height_direction == 1:
		label='$A_{z>0}/A_total$'
		areas_avg = [mean(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
		mytitle = 'Tile area (within $'+str(cutoff_distance)+'$ nm of protein) where $z>0$'
	else:
		label='$A_{z<0}/A_total$'
		areas_avg = [mean(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
		mytitle = 'Tile area (within $'+str(cutoff_distance)+'$ nm of protein) where $z<0$'
	areas_avg = [mean(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
	areas_err = [std(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
	ax.plot(cutoff_distance_sweep,areas_avg,color='b',lw=3,label=label)
	ax.errorbar(cutoff_distance_sweep,areas_avg,yerr=areas_err,color='b',fmt='-.',lw=3,label=label,capthick=4)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xlim((0,max(cutoff_distance_sweep)+1))
	ax.set_ylim((0.0,1))
	plt.xlabel('Distance from protein ($nm$)', labelpad = 10)
	ylabel1 =plt.ylabel('Area fraction', labelpad = 10)
	#plt.title(mytitle)
	ax.grid()
	plt.tight_layout()
	plt.savefig(pickles+'fig-'+sysname+'-tilefilter.area.sweep.png',dpi=500,
		bbox_extra_artists=[ylabel1])
	plt.close()
	
def batch_calculate_tilefilter_areas(make_figs=None,end=None,start=None,skip=None,framecount=None):
	'''Perform discretized area analysis for target frames.'''
	area_counts = []
	#---Make a directory for figures
	if make_figs:
		if not os.path.isdir(pickles+'/figs-'+sysname+'-tilefilter.snapshots'):
			os.mkdir(pickles+'/figs-'+sysname+'-tilefilter.snapshots')
	#---Loop over frames, and calculate areas and print figures
	nframes = len(mset.surf)
	if framecount == None:
		if end == None: end = nframes
		if start == None: start = 0
		if skip == None: skip = 1
	else:
		start = 0
		end = nframes
		skip = int(float(nframes)/framecount)
		skip = 1 if skip < 1 else skip
	result_data = MembraneData('tilefilter_area_v1')
	for fr in range(start,end,skip):
		if make_figs:
			result = lateral_discretize(fr,result='all')
			area_counts.append(result[0])
			blockplot(result[1]+result[2],zvals=zvals,colorcodes=colorcodes,
				filename=(pickles+'/figs-'+sysname+'-tilefilter.snapshots/fig%05d.png'%fr),show=False)
		else:
			result = lateral_discretize(fr,result='area')
			area_counts.append(lateral_discretize(fr))
			result_data.add(area_counts[-1],[fr])
	#---previously returned array(area_counts), but modified to use membraindata object
	return result_data

#---MAIN
#-------------------------------------------------------------------------------------------------------------

for ad in analysis_descriptors[analysis_plan]:
	(startpickle,protein_subset_slice,protein_pickle,suffix) = ad
	sysname = startpickle[24:-4]
	#---note: mset is global here
	mset = unpickle(pickles+startpickle)
	if protein_pickle != None:
		mset_protein = unpickle(pickles+protein_pickle)
		proteins_all = array(mset_protein.protein)
	else:
		proteins_all = array(mset.protein)
	proteins = proteins_all[:,protein_subset_slice]
	vecs = mean(mset.vecs,axis=0)
	cutoff = cutoff_distance*10/(vecs[0]/mset.griddims[0])
	print 'loaded '+startpickle
	print 'frame count = '+str(len(mset.surf[0]))
	result_data = batch_calculate_tilefilter_areas()
	result_data.addnote(['area_per_tile',product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))])
	pickle.dump(result_data,open(pickles+'pkl.tilefilter-areas.'+sysname+suffix+'.pkl','w'))


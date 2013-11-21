#!/usr/bin/python

#clrs = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors

#pts = array(A)
#pts = array(where(pts!=0)).T

#---mpl patch example
if 0:
	fig = plt.figure()
	ax = fig.add_subplot(111,aspect='equal')
	ax.set_xlim((0,33))
	ax.set_ylim((0,33))
	colors = 100*np.random.random(len(pts))
	patches = []
	for p in pts:
		rect = mpl.patches.Rectangle(tuple(p),1,1)
		patches.append(rect)
	p = mpl.collections.PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
	p.set_array(colors)
	ax.add_collection(p)
	plt.colorbar(p)
	plt.show()


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
		p = mpl.collections.PatchCollection(patches,cmap=brewer2mpl.get_map('Paired','qualitative',4,'reverse').mpl_colormap,alpha=1.0,edgecolors='none')
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

######## SIMPLEST METHOD

if 0:
	from membrainrunner import *;location='dark';execfile('locations-rpb.py')
	mset=unpickle(pickles+'pkl.avgstruct.membrane-v599.relevant.pkl')

if 0:
	#surf_discrete = [[(1 if mset.surf[10][i][j] > 0 else -1) for j in range(mset.griddims[1])] for i in range(mset.griddims[0])]
	mset.calculate_average_surface()
	surf_discrete = [[(1 if mset.surf_mean[i][j] > 0 else -1) for j in range(mset.griddims[1])] for i in range(mset.griddims[0])]
	#imshow(array(surf_discrete).T,origin='lower',interpolation='none');plt.show()
	#meshplot(mset.surf[10])
	#exeimshow(mset.surf[10],origin='lower');plt.show()

if 0:
	#---just plot protein and surface in imshow
	imshow(array(array(mset.rezipgrid(mset.protein[0]))!=0.0).T,origin='lower',interpolation='none',alpha=0.5)
	imshow(array(surf_discrete).T,origin='lower',interpolation='none',alpha=0.5)
	plt.show()
	
if 0:
	vecs = mean(mset.vecs,axis=0)
	cutoff = 50./(vecs[0]/mset.griddims[0])
	protpos = array(where(array(mset.rezipgrid(mset.protein[0]))!=0.0)).T
	#gridpos = [[i*vecs[0]/mset.griddims[0],j*vecs[1]/mset.griddims[1]] for i in range(mset.griddims[0]) for j in range(mset.griddims[1])]
	gridpos = [[i,j] for i in range(mset.griddims[0]) for j in range(mset.griddims[1])]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	#buf = [gridpos[j] for j in list(set([k for k in list(where(distmat[i]<cutoff)[0]) for i in range(len(distmat))]))]
	buf = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) for i in range(len(distmat))] for w in v]))]
	buf2 = array(scipy.sparse.coo_matrix(([1 for i in range(len(buf))],array(buf).T),dtype=int8,shape=(mset.griddims[0],mset.griddims[1])).todense())
	#array(A.todense())
	#imshow(array(surf_discrete).T,origin='lower',interpolation='none',alpha=0.2)
	#imshow(array(A).T,origin='lower',interpolation='none',alpha=0.2)
	#plt.show()
	
if 0:
	vecs = mean(mset.vecs,axis=0)
	cutoff = 50./(vecs[0]/mset.griddims[0])
	area_counts = []
	for fr in range(len(mset.surf)):
		print 'area check '+str(fr)
		surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] for i in range(mset.griddims[0]-1)])
		protpos = array(where(array(mset.rezipgrid(mset.protein[fr]))!=0.0)).T
		gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
		distmat = scipy.spatial.distance.cdist(protpos,gridpos)
		bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) for i in range(len(distmat))] for w in v]))]
		buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
		area_counts.append([float(sum(surf_discrete==1)),float(sum(surf_discrete==-1)),float(sum(surf_discrete+buf==2)),float(sum(surf_discrete+buf==1)),float(sum(buf==1))])
	area_counts = array(area_counts)
	
if 0:
	fig, ax = plt.subplots(1)
	t = range(len(area_counts))
	posarea = array([area_counts[i,2] for i in range(len(area_counts))])
	negarea = array([area_counts[i,3] for i in range(len(area_counts))])
	ax.plot(negarea,'b-')
	ax.plot(posarea,'r-')
	ax.fill_between(t, negarea,posarea, facecolor='blue', alpha=0.5,where=negarea>posarea)
	ax.grid()
	plt.show()

if 0:
	zvals={-1:0,1:1,0:2,2:3}
	which_brewer_colors = [1,5,0,4]
	colorcodes = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]

	fr=1000
	surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] for i in range(mset.griddims[0]-1)])
	protpos = array(where(array(mset.rezipgrid(mset.protein[fr]))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) for i in range(len(distmat))] for w in v]))]
	buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())	
	blockplot(surf_discrete+buf,zvals=zvals,colorcodes=colorcodes)
	
if 0:
	'''Diagnostic plot of frame filter figure, mean configuration.'''
	result = lateral_discretize(fr,result='grid')
	blockplot(result[0]+result[1],zvals=zvals,colorcodes=colorcodes)
	surf_discrete = array([[(1 if mset.surf_mean[i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] 
		for i in range(mset.griddims[0]-1)])
	protpos = array(where(array(mset.rezipgrid(mean(proteins,axis=0)))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
		dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())

if 0:
	'''Generate figures and print results.'''
	mpl.rcParams.update({'font.size': 14})
	area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((2,1),(0,0))
	t = range(len(area_counts))
	posarea = array([area_per_tile*area_counts[i,2] for i in range(len(area_counts))])
	negarea = array([area_per_tile*area_counts[i,3] for i in range(len(area_counts))])
	ax.plot(negarea,'b-',label='$z<0$',lw=2)
	ax.plot(posarea,'r-',label='$z>0$',lw=2)
	ax.fill_between(t, negarea,posarea, facecolor=colorcodes[zvals[-1]], alpha=0.5,where=negarea>posarea)
	ax.fill_between(t, posarea,negarea, facecolor=colorcodes[zvals[1]], alpha=0.5,where=negarea<posarea)
	ax.set_xlim((0,len(area_counts)))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Frame', labelpad = 10)
	ylabel1 =plt.ylabel('Area ($nm^2$)', labelpad = 10)
	plt.title('Tile area (within $'+str(cutoff_distance)+'$ nm of protein) by midplane height')
	plt.legend()
	ax.grid()
	ax2 = plt.subplot2grid((2,1),(1,0))
	histplot1 = ax2.hist(negarea,color=colorcodes[zvals[-1]],alpha=0.8,label='$z<0$')
	histplot2 = ax2.hist(posarea,color=colorcodes[zvals[1]],alpha=0.8,label='$z>0$')
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.spines['bottom'].set_position(('outward', 20))
	ax2.spines['left'].set_position(('outward', 30))
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ylabel2 = plt.ylabel('Frames', labelpad = 10)
	plt.xlabel('Area ($nm^2$)', labelpad = 10)
	plt.legend()
	ax2.grid()
	plt.savefig(pickles+'result.fig.tilefilter.areas.'+systemprefix+'.png',bbox_extra_artists=[ylabel1,ylabel2],bbox_inches='tight')
	plt.tight_layout()
	plt.close()
	fp = open(pickles+'result.txt.tilefilter.'+systemprefix+'.txt','w')
	fp.write('mean total area: '+str(product(vecs[0:2])/100.)+'nm2\n')
	fp.write('mean positive area: '+str(mean(area_per_tile*area_counts[:,0]))+' nm2\n')
	fp.write('mean negative area: '+str(mean(area_per_tile*area_counts[:,1]))+' nm2\n')
	fp.write('mean protein/positive area: '+str(area_per_tile*mean(area_counts[:,2]))+' nm2\n')
	fp.write('mean protein/negative area: '+str(area_per_tile*mean(area_counts[:,3]))+' nm2\n')
	fp.close()


if 0:
	filename=(pickles+systemprefix+'-tile-filter-snapshots/fig.mean.png')
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
	
if 0:
	mpl.rcParams.update({'font.size': 14})
	area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((1,1),(0,0))
	#t = range(len(area_counts_sweep[0]))
	#posarea = array([area_per_tile*area_counts_sweep[0][i,2] for i in range(len(area_counts_sweep[0]))])
	#negarea = array([area_per_tile*area_counts_sweep[0][i,2] for i in range(len(area_counts_sweep[0]))])
	#ax.plot(posarea,'r-',label='$z>0$',lw=3)
	#ax.plot(negarea,'b-',label='$z<0$',lw=3)
	#ax.fill_between(t, negarea,posarea, facecolor=colorcodes[zvals[-1]], alpha=0.5,where=negarea>posarea)
	#ax.fill_between(t, posarea,negarea, facecolor=colorcodes[zvals[1]], alpha=0.5,where=negarea<posarea)
	#ax.set_xlim((0,len(area_counts_sweep[0])))
	if height_direction == 1:
		label='$A_{z>0}/A_total$'
		areas_avg = [mean(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
	else:
		label='$A_{z<0}/A_total$'
		areas_avg = [mean(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
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
	plt.xlabel('Distance from protein ($nm$)', labelpad = 10)
	ylabel1 =plt.ylabel('Area fraction', labelpad = 10)
	plt.title('Tile area (within $'+str(cutoff_distance)+'$ nm of protein) by midplane height')
	ax.grid()
	plt.tight_layout()

def view_figures_sweep(area_counts_sweep):	
	mpl.rcParams.update({'font.size': 14})
	area_per_tile = product(vecs[0:2])/100./((mset.griddims[0]-1)*(mset.griddims[1]-1))
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((1,1),(0,0))
	if height_direction == 1:
		label='$A_{z>0}/A_total$'
		areas_avg = [mean(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,2])/mean(i[:,4]) for i in area_counts_sweep]
	else:
		label='$A_{z<0}/A_total$'
		areas_avg = [mean(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
		areas_err = [std(i[:,3])/mean(i[:,4]) for i in area_counts_sweep]
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
	ax.set_ylim((0.,1))
	plt.xlabel('Distance from protein ($nm$)', labelpad = 10)
	ylabel1 =plt.ylabel('Area fraction', labelpad = 10)
	plt.title('Tile area (within $'+str(cutoff_distance)+'$ nm of protein) where $z>0$')
	ax.grid()
	plt.tight_layout()
	plt.savefig(pickles+'result.fig.tilefilter.area-sweep.'+systemprefix+'.png',dpi=300,
		bbox_extra_artists=[ylabel1])
	plt.close()

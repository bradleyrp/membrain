#!/usr/bin/python -i

if 0:
	#---memory intensive
	from membrainrunner import *

	location = ''
	execfile('locations.py')
	execfile('plotter.py')

	#---selections
	director_cgmd = ['name PO4','name C4A','name C4B']
	selector_cgmd = 'name PO4'
	cgmd_protein = 'name BB'
	sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']

	#---analysis plan
	analysis_plan = slice(None,None)
	analysis_descriptors = [
		('tmpdel/pkl.tilt.v701.md.part0003.60000-160000-200.pkl',
		'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
		'v701.md.part0003.60000-160000-200',
		'membrane-v701',director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
		('tmpdel/pkl.tilt.v700.md.part0002.100000-200000-200.pkl',
		'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',
		'v700.md.part0002.100000-200000-200',
		'membrane-v700',director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
		('tmpdel/pkl.tilt.v550.md.part0006.200000-300000-200.pkl',
		'pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',
		'v550.part0006.200000-300000-200',
		'membrane-v550',director_cgmd,selector_cgmd,None,slice(-1,None))]

	#---extra imports
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

	#---FUNCTIONS
	#-------------------------------------------------------------------------------------------------------------

	vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
	find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices 
		if x in simplex for indx in simplex if indx !=x))

	def curvcalc(z,lenscale):
		'''Calculate mean and Gaussian curvature directly.'''
		zy, zx  = numpy.gradient(z,lenscale)
		zxy, zxx = numpy.gradient(zx,lenscale)
		zyy, _ = numpy.gradient(zy,lenscale)
		H = (zx**2 + 1)*zyy - 2*zx*zy*zxy + (zy**2 + 1)*zxx
		H = -H/(2*(zx**2 + zy**2 + 1)**(1.5))
		K = ((zxx*zyy)-(zxy)**2)
		K = -K/(2*(zx**2 + zy**2 + 1)**(1.5))
		return [H,K]

	#---MAIN
	#-------------------------------------------------------------------------------------------------------------

	allcurvtilts = []
	allcurvms = []

	for ad in analysis_descriptors[analysis_plan]: 
		starttime = time.time()
		print 'Starting analysis job.'
		(tiltpickle,structpickle,sysname,tests,director,selector,protein_selection,trajno) = ad
		#---load
		mset_tilt = unpickle(pickles+tiltpickle)
		mset_struct = unpickle(pickles+structpickle)
		#---load original for the positions
		traj = trajectories[systems.index(tests)][trajno][0]
		mset = MembraneSet()
		gro = structures[systems.index(tests)]
		basename = traj.split('/')[-1][:-4]
		sel_surfacer = sel_aamd_surfacer
		print 'Accessing '+basename+'.'
		mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
			resolution='cgmd')
		mset.identify_monolayers(director,startframeno=0)
		#---calculate mean curvature
		mset_struct.calculate_average_surface()
		lenscale = mean(mset_struct.vecs,axis=0)[0]/10./(shape(mset_struct.surf[0])[0]+1)
		cdat = lenscale*curvcalc(mset_struct.surf[0],lenscale)[0]
		curvsm = []
		for i in range(len(mset_struct.surf)):
			curvsm.append(curvcalc(list(array(mset_struct.surf[i]).T),lenscale)[0])
		#---calculate Gaussian curvature	
		lenscale = mean(mset_struct.vecs,axis=0)[0]/10./(shape(mset_struct.surf[0])[0]+1)
		curvsk = []
		for i in range(len(mset_struct.surf)):
			curvsk.append(curvcalc(list(array(mset_struct.surf[i])),lenscale)[1])
		#---calculate tilt angles
		ngridpts = shape(curvsm)[1]
		taila = mset_tilt.getdata('tilts').get(['monolayer',0,'tail',0])
		tailb = mset_tilt.getdata('tilts').get(['monolayer',0,'tail',1])
		framemap = [i[0] for i in mset_tilt.store[0].label]
		tbms = []
		for fr in range(len(framemap)):
			print fr
			frameno = framemap[fr]
			if frameno in mset_struct.vecs_index:
				lenscale = mset_struct.vec(frameno)[0]/10./(shape(mset_struct.surf[frameno])[0]-1)
				tilts_binned = [[[] for i in range(ngridpts)] for j in range(ngridpts)]
				for l in range(len(mset.monolayer_residues[0])):
					pos = mset.universe.residues[mset.monolayer_residues[0][l]].selectAtoms(
						selector).coordinates()[0]
					tilts_binned[int(round(pos[0]/lenscale/10.))][int(round(pos[1]/lenscale/10.))].\
						append(mean([taila[fr][l],tailb[fr][l]],axis=0))
				tilts_binned_mean = [[mean(tilts_binned[i][j]) if tilts_binned[i][j] != [] else 0. 
					for i in range(ngridpts)] for j in range(ngridpts)]
				tbms.append(tilts_binned_mean)
		curvtilts = [[180./pi*tbms[i][j][k],curvsm[i][j][k]] for i in range(shape(tbms)[0])
			for j in range(shape(tbms)[1]) for k in range(shape(tbms)[2])]
		allcurvtilts.append(curvtilts)
		allcurvms.append(curvsm)
		del mset
		del mset_struct
		del mset_tilt

if 0:	
	
	from membrainrunner import *
	
	location = ''
	execfile('locations.py')
	execfile('plotter.py')

	allcurvtilts = array(pickle.load(open(pickles+'pkl.tilt.allcurvtilts.pkl','r')))
if 0:	
	xspan = (min([min(allcurvtilts[i][:,0]) for i in range(3)]),max([max(allcurvtilts[i][:,0]) for i in range(3)]))
	yspan = (min([min(allcurvtilts[i][:,1]) for i in range(3)]),max([max(allcurvtilts[i][:,1]) for i in range(3)]))
	yspan = (-1*max([-1*yspan[0],yspan[1]]),max([-1*yspan[0],yspan[1]]))
	
	xspan = (140,165)
	yspan = (-0.6,0.6)
	
	nbins = 50

	extrema = 500

	names = ['antiparallel','parallel','control']

	fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(12,6))	
	for panel in range(3):
		ticknums = 5
		roundlev = -1
		curvtilts = allcurvtilts[panel]
		#xspan = (min(array(curvtilts)[:,0]),max(array(curvtilts)[:,0]))
		#yspan = (min(array(curvtilts)[:,1]),max(array(curvtilts)[:,1]))
		print 'computing H'
		H, xedges, yedges = histogram2d(array(curvtilts)[:,0],array(curvtilts)[:,1],
			range=(xspan,yspan),bins=nbins)
		midx = (xedges[1:]+xedges[:-1])/2
		midy = (yedges[1:]+yedges[:-1])/2
		extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
		cmap = mpl.cm.jet
		#cmap = mpl.cm.binary
		cmap.set_bad(cmap(0),1.)
		axes[0][panel].imshow(H.T, extent=None, interpolation='nearest',aspect='equal',
			origin='lower',norm=None,cmap=cmap,vmax=extrema,vmin=0)
	
		xts = [int(i) for i in arange(round(xspan[0],-1),round(xspan[1],-1)+0.001,10)]
		xtsinds = argmin([[abs(i-j) for i in midx] for j in xts],axis=1)
		axes[0][panel].axes.set_xticks(xtsinds)
		axes[0][panel].axes.set_xticklabels(xts)
	
		yts = [round(i,2) for i in arange(round(yspan[0],2),round(yspan[1],2),0.1)]
		ytsinds = argmin([[abs(i-j) for i in midy] for j in yts],axis=1)
		axes[0][panel].axes.set_yticks(ytsinds)
		axes[0][panel].axes.set_yticklabels(yts)
		axes[0][panel].axes.set_title(names[panel])
		#if panel == 0:
		axes[0][panel].axes.set_ylabel(r'$\mathsf{H(nm^{-1})}$')
		axes[0][panel].axes.set_xlabel('tilt angle (degrees)')
	controlH = H
	
	extrema = 65
	
	for panel in range(2):
		ticknums = 5
		roundlev = -1
		curvtilts = allcurvtilts[panel]
		#xspan = (min(array(curvtilts)[:,0]),max(array(curvtilts)[:,0]))
		#yspan = (min(array(curvtilts)[:,1]),max(array(curvtilts)[:,1]))
		print 'computing H'
		H, xedges, yedges = histogram2d(array(curvtilts)[:,0],array(curvtilts)[:,1],range=(xspan,yspan),
			bins=nbins)
		midx = (xedges[1:]+xedges[:-1])/2
		midy = (yedges[1:]+yedges[:-1])/2
		extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
		cmap = mpl.cm.jet
		#cmap = mpl.cm.binary
		cmap = 'bwr'
		#cmap.set_bad(cmap(0),1.)
		axes[1][panel].imshow((H-controlH).T, extent=None, interpolation='nearest',aspect='equal',
			origin='lower',norm=None,cmap=cmap,vmax=extrema,vmin=-extrema)
	
		xts = [int(i) for i in arange(round(xspan[0],-1),round(xspan[1],-1)+0.001,10)]
		xtsinds = argmin([[abs(i-j) for i in midx] for j in xts],axis=1)
		axes[1][panel].axes.set_xticks(xtsinds)
		axes[1][panel].axes.set_xticklabels(xts)
	
		yts = [round(i,2) for i in arange(round(yspan[0],2),round(yspan[1],2),0.1)]
		ytsinds = argmin([[abs(i-j) for i in midy] for j in yts],axis=1)
		axes[1][panel].axes.set_yticks(ytsinds)
		axes[1][panel].axes.set_yticklabels(yts)
		#if panel == 0:
		axes[1][panel].axes.set_ylabel(r'$\mathsf{H(nm^{-1})}$')
		axes[1][panel].axes.set_xlabel('tilt angle (degrees)')

	#axes[1][2].imshow([[0 for i in range(shape(H)[0])] for j in range(shape(H)[0])],cmap=mpl.cm.binary)
	axes[1][2].axis('off')
	if 0:
		minval=-0.2
		maxval=0.2
	
		array(allcurvtilts[0][:,0])
		hist0,binedge0 = numpy.histogram(array(allcurvtilts[0][:,1]),bins=nbins,normed=True,range=(minval,maxval))
		hist1,binedge1 = numpy.histogram(array(allcurvtilts[1][:,1]),bins=nbins,normed=True,range=(minval,maxval))
		hist2,binedge2 = numpy.histogram(array(allcurvtilts[2][:,1]),bins=nbins,normed=True,range=(minval,maxval))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		mid2 = (binedge2[1:]+binedge2[:-1])/2	
		axes[1][2].plot(mid0,hist0,'bo-',alpha=1.,lw=2,label=names[0])
		axes[1][2].plot(mid1,hist1,'co-',alpha=1.,lw=2,label=names[1])
		axes[1][2].plot(mid2,hist2,'ko-',alpha=1.,lw=2,label=names[2])
		axes[1][2].legend()
	plt.savefig(pickles+'fig-lipid-tilts-by-curvature-v700.v701.v550'+'.png',dpi=400)
	plt.show()

		
if 1:
	nbins = 30
	minval=-0.5
	maxval=0.5
	fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(8,6))	
	hist0,binedge0 = numpy.histogram(array(allcurvtilts[0][:,1]),bins=nbins,normed=True,range=(minval,maxval))
	hist1,binedge1 = numpy.histogram(array(allcurvtilts[1][:,1]),bins=nbins,normed=True,range=(minval,maxval))
	hist2,binedge2 = numpy.histogram(array(allcurvtilts[2][:,1]),bins=nbins,normed=True,range=(minval,maxval))
	mid0 = (binedge0[1:]+binedge0[:-1])/2
	mid1 = (binedge1[1:]+binedge1[:-1])/2
	mid2 = (binedge2[1:]+binedge2[:-1])/2	
	axes.plot(mid0,hist0,'bo-',alpha=1.,lw=2,label=names[0])
	axes.plot(mid1,hist1,'go-',alpha=1.,lw=2,label=names[1])
	axes.plot(mid2,hist2,'ko-',alpha=1.,lw=2,label=names[2])
	axes.legend()
	axes.axes.set_xlabel(r'$\mathsf{H(nm^{-1})}$')
	plt.savefig(pickles+'fig-curvature-distributions-v700.v701.v550'+'.png',dpi=400)
	plt.show()
		
		
		
		
		
		
		

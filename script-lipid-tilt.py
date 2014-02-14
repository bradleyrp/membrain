#!/usr/bin/python

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = None
framecount = None
location = ''
execfile('locations.py')

#---Selections
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v599-select'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v550'],director_cgmd,selector_cgmd,None,slice(-1,None))]

#---methods
do_lipid_tilt = True
do_lipid_tilt_vs_curv = True
	
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

#------------TEMPORARY
if 0:
	for ad in analysis_descriptors[analysis_plan]: 
		starttime = time.time()
		print 'Starting analysis job.'
		(startpickle,nprots,sysname,label) = ad
		#---load
		mset = unpickle(pickles+startpickle)
		#---calculate mean curvature
		mset.calculate_average_surface()
		lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
		cdat = lenscale*curvcalc(mset.surf[0],lenscale)[0]
		curvsm = []
		for i in range(len(mset.surf)):
			curvsm.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[0])
		#---calculate Gaussian curvature	
		lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
		curvsk = []
		for i in range(len(mset.surf)):
			curvsk.append(curvcalc(list(array(mset.surf[i])),lenscale)[1])
		#---plots
		extremum = max([np.max(curvsm),np.abs(np.min(curvsm))])
		extremum = abs(np.mean(curvsm))+2*np.std(curvsm)
		numgridpts = shape(curvsk)[1]
		fig = plt.figure(figsize=(12,3))
		ax0 = plt.subplot2grid((1,3), (0,0))
		ax0.set_xlim((0,numgridpts))
		ax0.set_ylim((0,numgridpts))
		ax0.set_title('mean curvature')
		img0 = ax0.imshow(array(np.mean(curvsm,axis=0)).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		if nprots > 0:
			protlen = int(shape(mset.protein[0])[0]/nprots)
			protlenshow = protlen
			for protpts in [mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)]:
				protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in protpts])
				hull = scipy.spatial.ConvexHull(protpts)
				ax0.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'k-',lw=0.6)
				shifthullx = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][0] 
					for i in range(len(hull.vertices))]
				shifthully = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][1] 
					for i in range(len(hull.vertices))]
				ax0.plot(shifthullx,shifthully,'k-',lw=0.6)		
		cax = inset_axes(ax0,
		     width="5%",
		     height="100%",
		     bbox_transform=ax0.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img0,cax=cax)
		cax.tick_params(labelsize=8) 
		cax.set_ylabel(r'$\mathsf{H(nm^{-1})}$',fontsize=10)
		extremum = max([np.max(curvsk),np.abs(np.min(curvsk))])
		extremum = abs(np.mean(curvsk))+2*np.std(curvsk)
		ax1 = plt.subplot2grid((1,3), (0,1))
		ax1.set_title('Gaussian curvature')
		ax1.set_xlim((0,numgridpts))
		ax1.set_ylim((0,numgridpts))
		img = ax1.imshow(array(np.mean(curvsk,axis=0)).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		if nprots > 0:
			protlen = int(shape(mset.protein[0])[0]/nprots)
			protlenshow = protlen
			for protpts in [mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)]:
				protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in protpts])
				hull = scipy.spatial.ConvexHull(protpts)
				ax1.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'k-',lw=0.6)
				shifthullx = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][0] 
					for i in range(len(hull.vertices))]
				shifthully = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][1] 
					for i in range(len(hull.vertices))]
				ax1.plot(shifthullx,shifthully,'k-',lw=0.6)					
		cax = inset_axes(ax1,
		     width="5%",
		     height="100%",
		     bbox_transform=ax1.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\mathsf{K(nm^{-2})}$',fontsize=10)
		extremum = max([np.max(mset.surf_mean),np.abs(np.min(mset.surf_mean))])
		ax2 = plt.subplot2grid((1,3), (0,2))
		ax2.set_title('structure')
		ax2.set_xlim((0,numgridpts))
		ax2.set_ylim((0,numgridpts))
		img2 = ax2.imshow(array(array(mset.surf_mean).T).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		if nprots > 0:
			protlen = int(shape(mset.protein[0])[0]/nprots)
			protlenshow = protlen
			for protpts in [mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow] 
				for i in range(nprots)]:
				protpts = array([[i[0]*numgridpts/vecs[0],i[1]*numgridpts/vecs[1]] for i in protpts])
				hull = scipy.spatial.ConvexHull(protpts)
				ax2.plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'k-',lw=0.6)
				shifthullx = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][0] 
					for i in range(len(hull.vertices))]
				shifthully = [protpts[hull.vertices[(i+1)%len(hull.vertices)]][1] 
					for i in range(len(hull.vertices))]
				ax2.plot(shifthullx,shifthully,'k-',lw=0.6)		
		cax = inset_axes(ax2,
		     width="5%",
		     height="100%",
		     bbox_transform=ax2.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img2,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\mathsf{z(nm)}$',fontsize=10)
		plt.tight_layout()
		plt.savefig(pickles+'fig-'+sysname+'-curvatures.mean.png',dpi=500,bbox_inches='tight')
		#plt.show()
		plt.cla()
		plt.close()

if do_lipid_tilt:
	starttime = time.time()
	print 'Starting analysis job.'
	for ad in analysis_descriptors[analysis_plan]:
		(tests,director,selector,protein_selection,trajno) = ad
		sysname = tests[0]
		residues = selector
		for t in range(len(tests)):
			testno = t
			print 'Running calculation: monolayer unstructured triangulation '+tests[t]+'.'
			for traj in trajectories[systems.index(tests[t])][trajno]:
				mset = MembraneSet()
				gro = structures[systems.index(tests[testno])]
				basename = traj.split('/')[-1][:-4]
				sel_surfacer = sel_aamd_surfacer
				print 'Accessing '+basename+'.'
				mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
					resolution='cgmd')
				mset.identify_monolayers(director,startframeno=0)
				#---frame selection header
				end = None
				start = None
				if framecount == None:
					if end == None: end = mset.nframes
					if start == None: start = 0
					if skip == None: skip = 1
				else:
					start = 0
					end = mset.nframes
					skip = int(float(mset.nframes)/framecount)
					skip = 1 if skip < 1 else skip
				#---initialize data objects
				result_data_surfnorms = MembraneData('surfnorms')
				result_data_tilts = MembraneData('tilts')
				if protein_selection != None:
					result_data_protdists = MembraneData('protdists')
				for fr in range(start,end,skip):
					print 'Calculating surface normals and distance to protein for frame: '+str(fr)+'.'
					mset.gotoframe(fr)
					#---get points
					topxyz = array([mean(mset.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
						for i in mset.monolayer_residues[0]])
					pts = mset.wrappbc(topxyz,mset.vec(fr),mode='grow')
					dtri = scipy.spatial.Delaunay(pts[:,0:2])
					point_permute = [[(i+j)%3 for i in range(3)] for j in range(3)]
					#---calculate triangle face normals
					trifaces = [[np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]) 
						for i in point_permute] for j in dtri.simplices]
					#---calculate simplex areas
					simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
					ptsareas = np.zeros(len(dtri.points))
					ptsnorms = np.zeros([len(dtri.points),3])
					print 'Summing points normals for each face.'
					#---calculate surface normals
					for s in range(len(dtri.simplices)):
						for p in range(3):
							ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
								if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
						ptsareas[dtri.simplices[s][p]] += simp_areas[s]
					ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
					result_data_surfnorms.add([[ptsnorms,[]],[ptsareas,[]]],[fr])
					#---calculate minimum distance to the protein
					if protein_selection != None:
						locations_protein = mset.universe.selectAtoms(protein_selection).coordinates()
						dmat = scipy.spatial.distance.cdist(locations_protein,topxyz)
						minprotdists = [min(i) for i in dmat.T]
						result_data_protdists.add([[minprotdists,[]]],[fr])
					#---calculate lipid tilts for each tail
					#---note that this uses the second and third items in the director as the tail atoms
					vecslipidsa = [mset.universe.residues[mset.monolayer_residues[0][i]].\
						selectAtoms(director_cgmd[1]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
					vecslipidsb = [mset.universe.residues[mset.monolayer_residues[0][i]].\
						selectAtoms(director_cgmd[2]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
					tiltsa = [arccos(np.dot(vecslipidsa[i],ptsnorms[i])/np.linalg.norm(
						vecslipidsa[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
					tiltsb = [arccos(np.dot(vecslipidsb[i],ptsnorms[i])/np.linalg.norm(
						vecslipidsb[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
					result_data_tilts.add([[tiltsa,[]],[tiltsb,[]]],[fr])
				mset.store.append(result_data_surfnorms)
				mset.store.append(result_data_tilts)
				if protein_selection != None:
					mset.store.append(result_data_protdists)			
				pickledump(mset,'pkl.tilt.'+sysname[9:14]+'.'+basename+'.pkl',directory=pickles)
				if erase_when_finished:
					del mset
	print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


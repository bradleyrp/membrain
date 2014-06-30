#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

if 'batch_override' not in globals():
	execfile('header-cgmd.py')
	execfile('header-meso.py')
	#---standard selection
	analysis_names = [
		'v530-40000-90000-50',
		'v531-40000-90000-50',
		'v532-40000-90000-50',
		'v509-40000-90000-50',
		'v510-40000-90000-50',
		'v511-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		'v514-22000-32000-10',
		'v515-20000-30000-10',
		'v614-120000-220000-200',
		'v2013-C_0-0.0369',
		][-1:]
	routine = [
		'calcs',
		]

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

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'calcs' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]

		#---prepare to process either an MD simulation or a mesoscale one
		if (analysis_descriptors[aname])['simtype'] == 'meso':
			mset = MembraneSet()
			mset.load_points_vtu(locate,extra_props='induced_cur',
				start=start,end=end,nbase=nbase,prefix='EQUIB-conf-')
			whichframes = None
			whichframes = slice(0,10)
		else:
			grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
			#---purposefully disabled looping over trajectory files so only one is assumed
			traj = trajfile[0]
			#---note: not checking to see if the pkl exists from this space
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
	
		result_data_points = MembraneData('cells')
		result_data_vmap = MembraneData('cells')
		result_data_areas = MembraneData('cells')
		result_data_surfnorms = MembraneData('surfnorms')
		if whichframes == None: whichframes = slice(None,None)
		#---loop over frames
		for frameno in range(mset.nframes)[whichframes]:
			status('status: processing frame ',i=frameno,looplen=mset.nframes)
			'''
			Todo for this section:
			calculate surface normals
			calculate lipid tilts
			calculate link lengths and identities
			does this require delaunay and voronoi?
			'''

			if (analysis_descriptors[aname])['simtype'] != 'meso':
				#---select the frame
				mset.gotoframe(frameno)
				if type(selector) == str:
					selector = [selector for i in range(len(mset.resnames))]
				#---collect the lipid positions
				#---needs meso here
				top = []
				for t in range(len(mset.resnames)):
					top.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid_abs[0][t]])
				bot = []
				for t in range(len(mset.resnames)):
					bot.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid[1][t]])
				pointgroups = [top,bot]
			else:
				mids = mset.xyzs[frameno]
				pointgroups = [mids]
				
			#---data collectors
			results_points = []
			results_vmap = []
			results_areas = []
			ptsnorms_monos = []
			areas_monos = []
			#---loop over monolayers
			for points in pointgroups:

				#---unwrap points and account for PBCs
				points_pbc = mset.wrappbc(points,vecs=mset.vec(frameno),mode='grow')

				#---AREAS
				#---compute voronoi mesh
				vmap = scipy.spatial.Voronoi(points_pbc[:,0:2])
				#---collect areas
				areas = []
				for p in range(len(points)):
					vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
					pairs = zip(vertices, vertices[1:]+vertices[0:1])
					areas.append(abs(sum(x1*y2-y1*x2 for (x1, y1), (x2, y2) in pairs)/2))
				#---storage is broken into parts for cell objects
				#---separately pickle the point positions, voronoi object, and areas
				results_points.append([points,[],[]])
				results_vmap.append([[],vmap,[]])
				results_areas.append([[],[],areas])

				#---NORMALS
				#---compute Delaunay
				dtri = scipy.spatial.Delaunay(points_pbc[:,0:2])
				point_permute = [[(i+j)%3 for i in range(3)] for j in range(3)]
				#---calculate triangle face normals
				trifaces = [[np.cross(points_pbc[j[i[1]]]-points_pbc[j[i[0]]],
					points_pbc[j[i[2]]]-points_pbc[j[i[0]]]) 
					for i in point_permute] for j in dtri.simplices]
				#---calculate simplex areas
				simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
				ptsareas = np.zeros(len(dtri.points))
				ptsnorms = np.zeros([len(dtri.points),3])
				#---calculate surface normals
				for s in range(len(dtri.simplices)):
					for p in range(3):
						ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
							if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
					ptsareas[dtri.simplices[s][p]] += simp_areas[s]
				ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
				ptsnorms_monos.append(ptsnorms)
				areas_monos.append(ptsareas)

				#---PROTEIN DISTANCES
				if 0 and protein_select != None:
					locations_protein = mset.universe.selectAtoms(protein_select).coordinates()
					dmat = scipy.spatial.distance.cdist(locations_protein,topxyz)
					minprotdists = [min(i) for i in dmat.T]

				#---TILTS
				#---note that this uses the second and third items in the director as the tail atoms
				if 0:
					vecslipidsa = [mset.universe.residues[mset.monolayer_residues[0][i]].\
						selectAtoms(director[1]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
					vecslipidsb = [mset.universe.residues[mset.monolayer_residues[0][i]].\
						selectAtoms(director[2]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
					tiltsa = [arccos(np.dot(vecslipidsa[i],ptsnorms[i])/np.linalg.norm(
						vecslipidsa[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
					tiltsb = [arccos(np.dot(vecslipidsb[i],ptsnorms[i])/np.linalg.norm(
						vecslipidsb[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]

			#---needs all add statements at this level, not the next level of indentation
			
			#---storage statements
			if 0: result_data_tilts.add([[tiltsa,[]],[tiltsb,[]]],[frameno])
			if 0: result_data_protdists.add([[minprotdists,[]]],[frameno])
			result_data_surfnorms.add([ptsnorms_monos,areas_monos],[frameno])
			result_data_points.add(results_points,[frameno])
			result_data_vmap.add(results_vmap,[frameno])
			result_data_areas.add(results_areas,[frameno])
			
			#---find a place to store the delaunay for later discrete clustering stuff
			
			#---dump everything
			
			datalist = [
				result_data_surfnorms,
				result_data_points,
				result_data_vmap,
				result_data_areas]

			for pklobj in datalist:
				for i in analysis_descriptors[aname]:
					pklobj.addnote([i,(analysis_descriptors[aname])[i]])
			
			for pklobj in datalist:
				mset.store = [pklobj]

			#---needs sysname 
				
			pickledump(mset,'pkl.cells.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)
				
				


#---LOADS NEEDS EDITED
#-------------------------------------------------------------------------------------------------------------

span_calcs = [
	'plot_span',
	'plot_span_angle',
	'plot_span_angle_summary',
	]
apl_calcs = [
	'plot_apl',
	'plot_apl_all_lipids',
	]

if routine in apl_calcs and 'msets' not in globals():
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.cells.'+\
			'areas-'+'-'.join(flaglist)+('-' if len(flaglist) > 0 else '')+\
			specname_pickle(sysname,trajfile[0])+'.pkl')
		msets.append(mset)
		
if routine in span_calcs and 'msets' not in globals():

	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.headspan-headangle2.'+specname_pickle(sysname,trajfile[0])+'.pkl')
		msets.append(mset)


#---DEPRECATED CODE BELOW
#-------------------------------------------------------------------------------------------------------------

#---generate voronoi cells
if 'xxxvoronoi' in routine:
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director,startframeno=0)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			result_data = MembraneData('cells')
			if whichframes == None: whichframes = slice(None,None)
			#---calcalate areas per lipid
			for frameno in range(mset.nframes)[whichframes]:
				print frameno
				#---select the frame
				mset.gotoframe(frameno)
				if type(selector) == str:
					selector = [selector for i in range(len(mset.resnames))]
				#---collect the lipid positions
				top = []
				for t in range(len(mset.resnames)):
					top.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid[0][t]])
				bot = []
				for t in range(len(mset.resnames)):
					bot.extend([mean(mset.universe.residues[i].selectAtoms(selector[t]).coordinates(),
						axis=0) for i in mset.monolayer_by_resid[1][t]])
				#---calculate areas for each monolayer
				results = []
				discard = False
				for points in [top,bot]:
					points_pbc = mset.wrappbc(points,vecs=mset.vec(frameno),mode='grow')
					try: vmap = scipy.spatial.Voronoi(points_pbc[:,0:2])
					except: discard = True
					areas = []
					for p in range(len(points)):
						vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
						pairs = zip(vertices, vertices[1:]+vertices[0:1])
						areas.append(abs(sum(x1*y2-y1*x2 for (x1, y1), (x2, y2) in pairs)/2))
					results.append([points,vmap,areas])
					del points,vmap,areas
				if not discard: result_data.add(results,[frameno])
				#---add details to the store
				for i in analysis_descriptors[aname]:
					result_data.addnote([i,(analysis_descriptors[aname])[i]])
				result_data.addnote(['selector',selector])
				result_data.addnote(['director',director])
			mset.store.append(result_data)
			#---Save the data
			pickledump(mset,'pkl.cells.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)

#---calculate surface normals
if 'xxxsurfnorms' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director,startframeno=0)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			result_data = MembraneData('cells')
			if whichframes == None: whichframes = slice(None,None)
			mset.identify_monolayers(director,startframeno=0)
			#---initialize data objects
			result_data_surfnorms = MembraneData('surfnorms')
			result_data_tilts = MembraneData('tilts')
			if protein_select != None:
				result_data_protdists = MembraneData('protdists')
			for fr in range(mset.nframes)[whichframes]:
				if fr % 100 == 0: print 'status: frame = '+str(fr)
				mset.gotoframe(fr)
				ptsnorms_monos = []
				areas_monos = []
				for mono in range(2):
					#---get points
					xyz = array([mean(mset.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
						for i in mset.monolayer_residues[mono]])
					pts = mset.wrappbc(xyz,mset.vec(fr),mode='grow')
					dtri = scipy.spatial.Delaunay(pts[:,0:2])
					point_permute = [[(i+j)%3 for i in range(3)] for j in range(3)]
					#---calculate triangle face normals
					trifaces = [[np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]) 
						for i in point_permute] for j in dtri.simplices]
					#---calculate simplex areas
					simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
					ptsareas = np.zeros(len(dtri.points))
					ptsnorms = np.zeros([len(dtri.points),3])
					#---calculate surface normals
					for s in range(len(dtri.simplices)):
						for p in range(3):
							ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
								if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
						ptsareas[dtri.simplices[s][p]] += simp_areas[s]
					ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
					ptsnorms_monos.append(ptsnorms)
					areas_monos.append(ptsareas)
				result_data_surfnorms.add([ptsnorms_monos,areas_monos],[fr])
			for i in analysis_descriptors[aname]:
				result_data_surfnorms.addnote([i,(analysis_descriptors[aname])[i]])
			mset.store.append(result_data_surfnorms)
			#---Save the data
			pickledump(mset,'pkl.surfnorms.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)
			if erase_when_finished:
				del mset
				
#---calculate tilts
if 'xxxtilts' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director,startframeno=0)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			result_data = MembraneData('cells')
			if whichframes == None: whichframes = slice(None,None)
			mset.identify_monolayers(director,startframeno=0)
			#---initialize data objects
			result_data_surfnorms = MembraneData('surfnorms')
			result_data_tilts = MembraneData('tilts')
			if protein_select != None:
				result_data_protdists = MembraneData('protdists')
			for fr in range(mset.nframes)[whichframes]:
				if fr % 100 == 0: print 'status: frame = '+str(fr)
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
				#---calculate surface normals
				for s in range(len(dtri.simplices)):
					for p in range(3):
						ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
							if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
					ptsareas[dtri.simplices[s][p]] += simp_areas[s]
				ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
				result_data_surfnorms.add([[ptsnorms,[]],[ptsareas,[]]],[fr])
				#---calculate minimum distance to the protein
				if protein_select != None:
					locations_protein = mset.universe.selectAtoms(protein_select).coordinates()
					dmat = scipy.spatial.distance.cdist(locations_protein,topxyz)
					minprotdists = [min(i) for i in dmat.T]
					result_data_protdists.add([[minprotdists,[]]],[fr])
				#---calculate lipid tilts for each tail
				#---note that this uses the second and third items in the director as the tail atoms
				vecslipidsa = [mset.universe.residues[mset.monolayer_residues[0][i]].\
					selectAtoms(director[1]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
				vecslipidsb = [mset.universe.residues[mset.monolayer_residues[0][i]].\
					selectAtoms(director[2]).coordinates()[0]-topxyz[i] for i in range(len(topxyz))]
				tiltsa = [arccos(np.dot(vecslipidsa[i],ptsnorms[i])/np.linalg.norm(
					vecslipidsa[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
				tiltsb = [arccos(np.dot(vecslipidsb[i],ptsnorms[i])/np.linalg.norm(
					vecslipidsb[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
				result_data_tilts.add([[tiltsa,[]],[tiltsb,[]]],[fr])
			for i in analysis_descriptors[aname]:
				result_data_surfnorms.addnote([i,(analysis_descriptors[aname])[i]])
				result_data_tilts.addnote([i,(analysis_descriptors[aname])[i]])
				if protein_select != None: result_data_protdists.addnote([i,(analysis_descriptors[aname])[i]])
			if protein_select != None:
				mset.store.append(result_data_protdists)
			mset.store.append(result_data_surfnorms)
			mset.store.append(result_data_tilts)
			#---Save the data
			pickledump(mset,'pkl.tilts.'+specname_pickle(sysname,traj)+'.pkl',directory=pickles)
			if erase_when_finished:
				del mset

#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = None
framecount = 3
location = ''
execfile('locations.py')

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v599-select'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
	(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None))]
	
#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices 
	if x in simplex for indx in simplex if indx !=x))

#---MAIN
#-------------------------------------------------------------------------------------------------------------

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
		
				#toptailxyz = array([mean(mset.universe.residues[i].selectAtoms(''.join([i+' or ' 
				#	for i in director[1:-1]]+[director[-1]])).coordinates(),axis=0) 
				#	for i in mset.monolayer_residues[0]])
				#print 'calculating angle'
				#vecslipids = [toptailxyz[i]-topxyz[i] for i in range(len(topxyz))]
				#print 1./60.*(time.time()-starttime)
				#meshpoints(array(ptsnorms)+array(pts))
				#meshplot(pts)
				#---calculate tilt angles
				#angles = [1./pi*arccos(np.dot(vecslipids[i],ptsnorms[i])/np.linalg.norm(
				#	vecslipids[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
				#plotthat = [[topxyz[i][0],topxyz[i][1],50*angles[i]] for i in range(len(topxyz))]
				#meshplot(plotthat)
				#areas = [1./3.*sum() for i in range(len(topxyz))]
				#[1./pi*arccos(np.dot(vecslipids[i],ptsnorms[i])/np.linalg.norm(
				#	vecslipids[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
				# adding a section for calculating minimum distance to a protein here
				#print 1./60.*(time.time()-starttime)
				#locations_protein = mset.universe.selectAtoms(protein_selection).coordinates()
				#locations_lipids = topxyz
				#dmat = scipy.spatial.distance.cdist(locations_protein,locations_lipids,metric=torus_norm)
				#dmat = scipy.spatial.distance.cdist(locations_protein,locations_lipids)
				#minprotdists = [min(i) for i in dmat.T]
				#print 1./60.*(time.time()-starttime)
				#result_data.add([[angles,[]],[ptsareas,[]]],[fr])
				#result_data_position.add([[topxyz,[]],[minprotdists,[]]],[fr])
				#if protein_selection != None:
				#	mset.protein.append(mset.universe.selectAtoms(protein_selection).coordinates())
				#	mset.protein_index.append(fr)
			mset.store.append(result_data_surfnorms)
			mset.store.append(result_data_protdists)
			mset.store.append(result_data_tilts)
			pickledump(mset,'pkl.tilt.'+sysname[9:14]+'.'+basename+'.pkl',directory=pickles)
	
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


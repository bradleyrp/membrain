#!/usr/bin/python -i

if 0:

	from membrainrunner import *

	#---Settings
	#-------------------------------------------------------------------------------------------------------------

	#---Analysis parameters
	skip = None
	framecount = None
	location = ''
	execfile('locations.py')

	#---Selections
	director_symmetric = ['name P','name C218','name C318']
	director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
			'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
	selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

	#---Analysis plan
	analysis_descriptors = [
		(['membrane-v510'],'all',director_symmetric,-1,
			[['resname DOPC and name P','name CL','DOPC P-to-CL'],
			['resname DOPC and name P','name MG','DOPC P-to-MG']]),
		(['membrane-v530'],'all',director_asymmetric,-1,
			[['resname DOPS and name P','resname PI2P and name P','DOPS-PIP2']]),
		(['membrane-v509'],'all',director_symmetric,-1,
			[['resname DOPS and name P','resname PI2P and name P','DOPS-PIP2']])]

	#---Functions
	#-------------------------------------------------------------------------------------------------------------

	(tests,selector,director,trajno,pairs) = analysis_descriptors[-1]
	t = 0
	traj = [trajectories[systems.index(tests[t])][trajno]][-1]
	testno = 0
	pair = pairs[0]

#def analyze_ion_distributions(testno,traj):
	'''Compute the lipid-lipid g(r).'''
	mset = MembraneSet()
	#---Load the trajectory
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]
	sel_surfacer = sel_aamd_surfacer
	print 'Accessing '+basename+'.'
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
		resolution='aamd')
	mset.identify_monolayers(director)
	self = mset
	
if 0:
	print 'code here'
	
	for lnum in range(2):
		allselect_lipids = self.universe.selectAtoms(pair[lnum])
		for mononum in range(2):
			validresids = list(set.intersection(set(self.monolayer_residues[mononum]),set([i-1 for i in allselect_lipids.resids()])))
			self.selections.append(sum([allselect_lipids.residues[list(allselect_lipids.resids()).index(i+1)].selectAtoms(pair[lnum]) for i in validresids]))

vecs = mset.vec(0)
def torus_norm(x1,x2):
	'''Provides the norm on a torus, given its dimensions.'''
	temp = x1-x2
	#if vecs == 0: vecs = self.vecs[0]
	if len(shape(temp)) == 1:
		#---checked documentation and tested cdist. ds and rpb found that cdist uses this option
		#---below is the re-written version for periodic boundary conditions
		#return sqrt(((x1-x2)**2).sum(axis=0)) 
		stddist = x1-x2
		stddist = array([abs(stddist[i])-(0.5 if stddist[i] > 0.5*vecs[i] else 0)*vecs[i] for i in range(len(x1))])
		return sqrt((stddist**2).sum(axis=0))
	elif len(shape(temp)) == 2:
		filt = array([[min([abs(dimswitch[d]*(temp[d][j]-vecs[d]*i)) for i in [-1,0,1]]) \
			for j in range(len(temp[d]))] for d in range(len(temp))])
		return sqrt(((filt)**2).sum(axis=0)) 
	elif len(shape(temp)) == 3:
		filt = array([[[min([abs((temp[d][j][k]-vecs[d]*i)) for i in [-1,0,1]]) \
			for k in range(len(temp[d][j]))] for j in range(len(temp[d]))] for d in range(len(temp))])
		return sqrt(((filt)**2).sum(axis=0)) 
	else:
		print 'Error: wrong dimensionality of the inputs to torus_norm.'
		return 0
		
vecs = mset.vec(0)
def crap_norm(x1,x2):
	'''Provides the norm on a torus, given its dimensions.'''
	temp = x1-x2
	#if vecs == 0: vecs = self.vecs[0]
	if len(shape(temp)) == 1:
		return sqrt(((x1-x2)**2).sum(axis=0))
	elif len(shape(temp)) == 2:
		filt = array([[min([abs(dimswitch[d]*(temp[d][j]-vecs[d]*i)) for i in [-1,0,1]]) \
			for j in range(len(temp[d]))] for d in range(len(temp))])
		return sqrt(((filt)**2).sum(axis=0)) 
	elif len(shape(temp)) == 3:
		filt = array([[[min([abs((temp[d][j][k]-vecs[d]*i)) for i in [-1,0,1]]) \
			for k in range(len(temp[d][j]))] for j in range(len(temp[d]))] for d in range(len(temp))])
		return sqrt(((filt)**2).sum(axis=0)) 
	else:
		print 'Error: wrong dimensionality of the inputs to torus_norm.'
		return 0
	
if 1:

	frameno = 0

	pts1 = array(self.get_points(frameno,selection_index=0))[:,0:2]
	pts2 = array(self.get_points(frameno,selection_index=0))[:,0:2]

	dmat = scipy.spatial.distance.cdist(pts1,pts2,metric=torus_norm)
	dmat2 = scipy.spatial.distance.cdist(pts1,pts2)
	
	if 1:
		hist0 = numpy.histogram([i for j in dmat for i in j])
		#plt.bar(hist0[1][:-1],2*hist0[0],color='b',alpha=0.5)
		plt.plot(hist0[1][:-1],1*hist0[0],'o-')
		plt.show()

	
	
	'''
	#---Batches of g(r)-voronoi calculations
	for pair in pairs:
		print pair
		mset.batch_gr_lipid_lipid([pair[0],pair[1]],framecount=framecount,skip=skip,
			label=pair[2],mode='voronoi_bin',monolayer_rep='P')
	#---Save the data
	pickledump(mset,'pkl.gr-vornoi.'+tests[testno]+'.'+basename+'.pkl')
	return mset
	'''

#---MAIN
#-------------------------------------------------------------------------------------------------------------

'''
starttime = time.time()
print 'Starting analysis job.'
for ad in analysis_descriptors:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,selector,director,trajno,pairs) = ad
	for t in range(len(tests)):
		print 'Running calculation: monolayer unstructured triangulation '+tests[t]+'.'
		for traj in [trajectories[systems.index(tests[t])][trajno]]:
			#---Run the analysis function on the desired system
			mset = analyze_ion_distributions(t,traj)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'
'''

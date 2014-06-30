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
			[['resname PI2P and name P','resname PI2P and name P','DOPS-DOPS']])]

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
		stddist = array([abs(stddist[i])-(0.5*vecs[i] if stddist[i] > 0.5*vecs[i] else 0) for i in range(len(x1))])
		if 0:		
			if any((stddist)!=abs(x1-x2)):
				print 'inside torus_norm'
				print temp
				print shape(temp)
				print (stddist)==abs(x1-x2)
				print x1
				print x2
				print 'stddist = '+str(stddist)
				print 'length = '+str(sqrt((stddist**2).sum(axis=0)))
		if sqrt((stddist**2).sum(axis=0)) < stddist[0] or sqrt((stddist**2).sum(axis=0)) < stddist[1]:
			print 'you done fucked u[p'
			print stddist
			print sqrt((stddist**2).sum(axis=0))
			print temp
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
	
if 0:

	frameno = 0

	pts1 = array(self.get_points(frameno,selection_index=0))[:,0:2]
	pts2 = array(self.get_points(frameno,selection_index=0))[:,0:2]

	dmat = scipy.spatial.distance.cdist(pts1,pts2,metric=torus_norm)
	dmat2 = scipy.spatial.distance.cdist(pts1,pts2)
	
	if 1:
		hist0,binedge0 = numpy.histogram([i for j in dmat for i in j],range=(0,250),bins=200)
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		hist1,binedge1 = numpy.histogram([i for j in dmat2 for i in j],range=(0,250),bins=200)
		mid1 = (binedge1[1:]+binedge1[:-1])/2
		binwidth = mid0[1]-mid0[0]
		areas = [((mid0[i]+binwidth/2.)**2-(mid0[i]-binwidth/2.)**2)*pi for i in range(len(mid0))]
		
		plt.plot(mid0,hist0,'o-',c='r',label='torus')
		plt.plot(mid1,hist1,'o-',c='b',label='reg')
		plt.plot(mid0,areas,'o-',c='g')
		plt.xlim((0,300))
		plt.legend()
		
		plt.show()
		
if 0:	
	
	pts1 = array(self.get_points(frameno,selection_index=0))[:,0:2]
	pts2 = array(self.get_points(frameno,selection_index=0))[:,0:2]

	smalldists = []
	images1 = []
	for i in range(len(pts1)):
		for j in range(len(pts2)):
			v = array([pts2[i][0]-pts1[j][0],pts2[i][1]-pts1[j][1]])
			#---images1
			#if abs(v[0]) <= 0.5*vecs[0] and abs(v[1]) <= 0.5*vecs[1]:
			smalldists.append(scipy.linalg.norm(v))
			for shift in [[-1,-1],[-1,0],[0,-1],[0,1],[1,0],[1,1],[1,-1],[-1,1]]:
				images1.append(scipy.linalg.norm(v+shift*vecs[0:2]))
	
	
if 1:
	allcurvs = []
	#frameno = 0
	for frameno in range(0,300,10):
		pts1 = array(self.get_points(frameno,selection_index=0))[:,0:2]
		pts2 = array(self.get_points(frameno,selection_index=0))[:,0:2]

		points = pts2
		dims=[0,1]
		ans = []
		for p in points:
			for tr in [[i,j] for i in arange(-2,2+1,1) for j in arange(-2,2+1,1) if not (i == 0 and j == 0)]:
				ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(2)])
		pts2pbc = concatenate((points,array(ans)))

		'''
		points = pts1
		dims=[0,1]
		ans = []
		for p in points:
			for tr in [[1,0],[0,1],[-1,0],[0,-1],[-1,1],[-1,-1],[1,1],[1,-1]]:
				ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(2)])
		pts1pbc = concatenate((points,array(ans)))
		points = pts2
		dims=[0,1]
		ans = []
		for p in points:
			for tr in [[1,0],[0,1],[-1,0],[0,-1],[-1,1],[-1,-1],[1,1],[1,-1]]:
				ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(2)])
		pts2pbc = concatenate((points,array(ans)))
		'''

		cutoff = 2*(vecs[0] if vecs[0]<vecs[1] else vecs[1])
	
		sysarea = pi*cutoff**2

		dmat2 = scipy.spatial.distance.cdist(pts1,pts2pbc)
		binsizeabs = 4
		hist,binedge = numpy.histogram(array(dmat2[0])[array(dmat2[0])!=0.],range=(0.,int(cutoff)),bins=int(cutoff/binsizeabs))
		mid = (binedge[1:]+binedge[:-1])/2
		areas = [pi*binwidth*mid[i]*2 for i in range(len(binedge)-1)]
		areas = [pi*((binedge[i+1])**2-(binedge[i])**2) for i in range(len(binedge)-1)]
	
		histcomb = []
		for r in range(len(dmat2)):
			row = dmat2[r]
			hist,binedge = numpy.histogram(array(row)[array(row)!=0.],range=(0.,int(cutoff)),bins=int(cutoff/binsizeabs))
			histcomb.append(hist)
			#if r%5 == 0:
			#	plt.plot(mid,hist/areas/(len(dmat2)/(vecs[0]*vecs[1])),'-',label='name',alpha=0.5)
		grcurve = sum(histcomb,axis=0)/float(len(histcomb))

		'''
		dmat2flat = array(dmat2).flatten()
		dmat2flat = dmat2flat[dmat2flat<=cutoff]
		dmat2flat = dmat2flat[dmat2flat!=0.]
		'''

		allcurvs.append(grcurve)

	mid = (binedge[1:]+binedge[:-1])/2
	binwidth = mid[1]-mid[0]
	plt.plot(mid,mean(allcurvs,axis=0)/areas/(len(dmat2)/(vecs[0]*vecs[1])),'o-',c='r',label='reg')
	plt.xlim((0,100))
	plt.ylim((0,2))
	#plt.legend()
	plt.show()
	

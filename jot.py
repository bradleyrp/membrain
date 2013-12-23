#!/usr/bin/python -i

if 0:
	from membrainrunner import *
	location = ''
	execfile('locations.py')
	mset = unpickle(pickles+'tilttest.pkl')

if 0:
	tmp = mset.store[0]
	dat = tmp.get(['monolayer',0,'type','angle'])

if 0:
	frame = 132
	self = mset
	mset.gotoframe(frame)
	topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])

if 0:
	hist, bin_edges = histogram([i for j in dat for i in j])
	mid = (bin_edges[1:] + bin_edges[:-1])/2
	plt.plot(mid, hist, 'o-',linewidth=2)
	plt.show()
	
if 0:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	d2 = scipy.spatial.distance.cdist(topxyz[:,0:2],array(mset.protein[1])[:,0:2])
	mindists = np.min(d2,axis=1)
	H, xedges, yedges = histogram2d(mindists,dat[0][0],bins=30)
	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	ax.set_xticklabels(xedges)
	ax.set_yticklabels(yedges)
	plt.imshow(H, extent=None, interpolation='nearest',aspect='equal',cmap='Reds')
	plt.show()
	
if 0:
	frame = 132
	fig = plt.figure()
	ax = fig.add_subplot(111)
	self = mset
	mset.gotoframe(frame)
	topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])
	d2 = scipy.spatial.distance.cdist(topxyz[:,0:2],array(mset.protein[1])[:,0:2])
	mindists = np.min(d2,axis=1)
	H, xedges, yedges = histogram2d(mindists,dat[0][0],bins=30)
	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	ax.set_xticklabels(xedges)
	ax.set_yticklabels(yedges)
	plt.imshow(H, extent=None, interpolation='nearest',aspect='equal',cmap='Reds')
	plt.show()

if 0:
	tmp2 = [[topxyz[i,0],topxyz[i,1],np.mean(dat[0],axis=0)[i]] for i in range(len(dat[0][0]))]
	meshplot(tmp2,show='surf')
	
if 0:
	tmp2 = [[topxyz[i,0],topxyz[i,1],np.mean(dat[0],axis=0)[i]] for i in range(len(dat[0][0]))]
	dat3 = mset.rezipgrid(array(tmp2),vecs=mset.vecs[0],grid=(20,20))
	plt.imshow(dat3, extent=None, interpolation='nearest',aspect='equal',cmap='Greys')
	plt.show()
	
if 0:
	tmp = mset.store[0]
	dat = tmp.get(['monolayer',0,'type','angle'])
	frame = 0
	self = mset
	mset.gotoframe(frame)
	topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	numgridpts = 8
	vecs = mset.vecs[0]
	dat4 = []
	for fr in range(len(dat)):
		print fr
		tmp2 = [[topxyz[i,0],topxyz[i,1],dat[fr][0][i]] for i in range(len(dat[fr][0]))]
		dat4.append(mset.rezipgrid(array(tmp2),vecs=mset.vecs[0],grid=(numgridpts,numgridpts)))
	plt.imshow(mean(dat4,axis=0), extent=None, interpolation='nearest',aspect='equal',cmap='Greys')
	for pt in mset.protein[0]:	
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='r',
			alpha=1.0)
		ax.add_patch(circ)
	plt.show()
if 0:
	tmp = mset.store[0]
	dat = tmp.get(['monolayer',0,'type','area'])
	frame = 0
	self = mset
	mset.gotoframe(frame)
	topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	numgridpts = 6
	vecs = mset.vecs[0]
	dat4 = []
	for fr in range(len(dat)):
		print fr
		tmp2 = [[topxyz[i,0],topxyz[i,1],dat[fr][0][i]] for i in range(6400)]
		dat4.append(mset.rezipgrid(array(tmp2),vecs=mset.vecs[0],grid=(numgridpts,numgridpts)))
	plt.imshow(mean(dat4,axis=0), extent=None, interpolation='nearest',aspect='equal',cmap='Greys')
	for pt in mset.protein[0]:	
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2*1.5/64*numgridpts,color='r',
			alpha=1.0)
		ax.add_patch(circ)
	plt.show()

##############################################

if 1:	
	from membrainrunner import *
	location = ''
	execfile('locations.py')
	mset = unpickle(pickles+'tilttest.pkl')
	

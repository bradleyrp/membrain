#!/usr/bin/python -i

if 0:

	from membrainrunner import *

	location = ''
	execfile('locations.py')
	execfile('plotter.py')

	mset = unpickle(pickles+'pkl.tilt.membrane-v599-select.md.part0003.select.nosol.pkl')

'''
...

'''

if 1:
	areas = mset.store[0].get(['monolayer',0,'type','area'])
	angles = mset.store[0].get(['monolayer',0,'type','angle'])
	positions = array(mset.store[1].data)[:,0]

	fig = plt.figure()
	ax = fig.add_subplot(111)

	numgridpts = 6
	vecs = mset.vecs[0]

	#---coallate data
	binned_data = []
	for fr in range(len(angles)):
		print fr
		combined_data = [[positions[fr][i,0],positions[fr][i,1],angles[fr][0][i]] for i in range(len(angles[0][0]))]
		binned_data.append(mset.rezipgrid(array(combined_data),vecs=mset.vecs[0],grid=(numgridpts,numgridpts)))
	plt.imshow(mean(binned_data,axis=0), extent=None, interpolation='nearest',aspect='equal',cmap='Greys')

	#---print one protein
	for pt in mset.protein[0]:
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2*1.5/64*numgridpts,color='r',
			alpha=1.0)
		ax.add_patch(circ)
	plt.show()



		

'''
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
'''

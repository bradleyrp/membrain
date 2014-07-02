#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

'''
'polyangles':{
	'struct':{'frame':0,'angles':1},
	'struct_opts':None,
	'description':'stores the angles along a polymer chain',
	},
'''

if 0:
	mset = MembraneSet()
	mset.load_trajectory((
		'/home/rpb/worker/demodat/system-input.gro',
		'/home/rpb/worker/demodat/md.part0005.skip10.xtc'
		))

	result_data = MembraneData('polyangles')
	for fr in range(mset.nframes):
		status('status: angle calculation, frame = ',i=fr,looplen=mset.nframes)
		mset.gotoframe(fr)
		#---alpha carbon positions
		pts = mset.universe.selectAtoms('name CA').coordinates()
		#---vectors between Ca along the protein backbone
		vecs = [pts[j+1]-pts[j] for j in range(0,len(pts)-1)]
		angle = [
			arccos(dot(vecs[i],vecs[i+1])/\
			linalg.norm(vecs[i])/\
			linalg.norm(vecs[i+1]
			)) for i in range(len(vecs)-1)]
		result_data.add(angle,[fr])
	mset.store.append(result_data)
	pickledump(mset,pickles+'mypkl2.pkl')
	del mset
	

if 1:
	mset = unpickle(pickles+'mypkl2.pkl')
	
	angles = mset.getdata('polyangles').data
	fig = plt.figure()
	for angle in angles:
		plt.plot(range(len(angle)),angle,c='k',alpha=0.2)
	plt.plot(range(len(angles[0])),mean(angles,axis=0),color='r')
	plt.show()
	


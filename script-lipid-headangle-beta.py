#!/usr/bin/python -i

if 1:

	from membrainrunner import *

	location = ''
	execfile('locations.py')

	# c2 P C14

	#---Analysis plan
	analysis_plan = slice(-1,None)
	analysis_descriptors = [
		(['membrane-v534'],None,'all',None,slice(-1,None),[60000,160000,200]),
		(['membrane-v533'],None,'all',None,slice(-1,None),[60000,160000,200])]

	t,testno = 0,0
	(tests,director,selector,protein_select,trajno,timeslice) = analysis_descriptors[-2]
	traj = trajectories[systems.index(tests[t])][trajno][-1]
	gro = structures[systems.index(tests[testno])]
	basename = traj.split('/')[-1][:-4]


	mset = MembraneSet()
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='cgmd')

if 1:
	angles = []
	for fr in range(100):
		print fr
		mset.gotoframe(fr)
		pts = mset.universe.selectAtoms('resname P35P and (name C2 or name P or name C14)').coordinates()
		angles.append([arccos(np.dot(pts[i]-pts[i+1],pts[i+2]-pts[i+1])/np.linalg.norm(pts[i]-pts[i+1])/np.linalg.norm(pts[i+2]-pts[i+1])) for i in range(0,120,3)])
	
if 1:
	hist0 = numpy.histogram([i for j in 180./pi*array(angles) for i in j],bins=40)
	plt.plot(hist0[1][:-1],1*hist0[0],'o-')
	plt.show()

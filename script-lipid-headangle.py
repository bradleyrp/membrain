#!/usr/bin/python -i

if 1:

	from membrainrunner import *

	location = ''
	execfile('locations.py')

	# c2 P C14

	#---Analysis plan
	analysis_plan = slice(-1,None)
		
	analysis_descriptors = {
	'v533-lipid-headangle':
		{'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-headangle',
		'trajsel':'md.part0012.xtc'},
	'v534-lipid-headangle':
		{'sysname':'membrane-v534',
		'sysname_lookup':'membrane-v534-headangle',
		'trajsel':'md.part0012.xtc'},
	}
	analysis_names = ['v533-lipid-headangle','v534-lipid-headangle']

	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')

		angles = []
		for fr in range(100):
			print fr
			mset.gotoframe(fr)
			pts = mset.universe.selectAtoms('resname P35P and (name C2 or name P or name C14)').coordinates()
			angles.append([arccos(np.dot(pts[i]-pts[i+1],pts[i+2]-pts[i+1])/np.linalg.norm(pts[i]-pts[i+1])/np.linalg.norm(pts[i+2]-pts[i+1])) for i in range(0,120,3)])
	

		hist0 = numpy.histogram([i for j in 180./pi*array(angles) for i in j],bins=40)
		plt.plot(hist0[1][:-1],1*hist0[0],'o-')
		plt.show()

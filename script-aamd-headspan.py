#!/usr/bin/python -i

from membrainrunner import *

from scipy import spatial

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
analysis_plan = slice(0,2)
analysis_descriptors = [
	(['membrane-v534'],['Cal'],['POPC','CHL1','DOPE','DOPS','P35P'],
		'all',director_asymmetric,slice(-1,None),
		'resname P35P and (name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'P35P'),
	(['membrane-v533'],['Mg'],['POPC','CHL1','DOPE','DOPS','P35P'],
		'all',director_asymmetric,slice(-1,None),
		'resname P35P and (name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'P35P'),
	(['membrane-v531-headspan'],['Mg'],['PI2P'],
		'all',director_asymmetric,slice(-1,None),
		'resname PI2P and (name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'PI2P'),
	(['membrane-v532-headspan'],['Cal'],['PI2P'],
		'all',director_asymmetric,slice(-1,None),
		'resname PI2P and (name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'PI2P')]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Compute P35P spans
tf_areas = []
starttime = time.time()
#---loop over analysis descriptors
for ad in analysis_descriptors[analysis_plan]:
	(tests,ionnames,residues,selector,director,trajno,selstring,targetresname) = ad
	#---loop over tests within the descriptor
	for testno in range(len(tests)):
		#---loop over specified trajectories
		for traj in trajectories[systems.index(tests[testno])][trajno]:
			mset = MembraneSet()
			gro = structures[systems.index(tests[testno])]
			basename = traj.split('/')[-1][:-4]
			print 'Accessing '+basename+'.'
			mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director,startframeno=0)
			mset.identify_residues(residues)
			result_data = MembraneData('cells')
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
			u_tf = mset.universe
			frameskip = 1
			tf_distance = []
			target_residues_abs = mset.monolayer_by_resid[0][mset.resnames.index(targetresname)]
			target_residues_abs = range(40)
			for frame in mset.universe.trajectory[::frameskip]:
				print frame
				selection = [mset.universe.residues[i].selectAtoms(selstring).coordinates() for i in target_residues_abs]
				tf_distance.append([max(spatial.distance.pdist(selection[i])) for i in range(40)])
			tf_areas.append([j*j for i in tf_distance for j in i])
			if erase_when_finished:
				del mset

pickle.dump(tf_areas,open(pickles+'pkl.headspan.v531.v532.part0010.25500-27500-4.pkl','w'))

if 0:
	fig = plt.figure()
	ax = plt.subplot(111)
	for sys in range(len(tf_areas)):
		hist0 = numpy.histogram(tf_areas[sys],bins=50)
		c = clrs[sys%len(clrs)]
		ax.plot(hist0[1][1:],hist0[0],'o-',lw=2,color=c)
	plt.savefig('/home/rpb/tmp.png')
	plt.show()


#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
analysis_descriptors = {
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl'}}
analysis_names = ['v511-30000-80000-100']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if mset.universe == []:
	print 'loading trajectory'
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset_surf = unpickle(pickles+structure_pkl)
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='cgmd')
		checktime()

#---get ion positions and times
if 0:
	print 'getting ions'
	clock = []
	ionspos = []
	ion_select = mset.universe.selectAtoms('name Cal')
	whichframes = range(len(mset.universe.trajectory))
	for fr in whichframes:
		print fr
		mset.gotoframe(fr)
		ionspos.append(ion_select.coordinates())
		clock.append(mset.universe.trajectory[fr].time)
	vecs=mean(mset_surf.vecs,axis=0)
	ionspos = array(ionspos)[:-1]
	ions_traj = []
	for ind in range(shape(ionspos)[1]):
		print ind
		course = array(ionspos)[:,ind]
		#---three-line handling PBCs
		hoplistp = (course[1:]-course[:-1])>array(mset_surf.vecs)[1:]/2.
		hoplistn = (course[:-1]-course[1:])>array(mset_surf.vecs)[1:]/2.
		course_nojump = course[1:]-(cumsum(1*hoplistp-1*hoplistn,axis=0))*array(mset_surf.vecs)[1:]
		ions_traj.append(course_nojump)
	ions_traj=array(ions_traj)
	nions = len(ions_traj)

#---GENERAL METHOD
if 0:
	center = mean(mset_surf.surf_position)
	zones = [[0,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]
	zones = zones + array(center)
	occupancy = 0.8
	alldists_filt_zone = []
	diffusion_direction = '2d'
	if diffusion_direction == '2d':
		dimslice = slice(0,2)
		msd_factor = 2
	elif diffusion_direction == 'z':
		dimslice = slice(2,3)
		msd_factor = 1
	else:
		dimslice = slice(None,None)
		msd_factor = 3
	for zone in zones:
		alldists_filt = []
		#---deprecated pre-filtering for zones, but mean z is almost always near bilayer
		if 0:
			ions_in_zone = where((1*(meanz>zone[0])+1*(meanz<zone[1]))==2)[0]
			print ions_in_zone
		else:
			ions_in_zone = range(nions)
		#---locate whether ions are in the zone at each time step
		for ind in ions_in_zone:
			print 'ion = '+str(ind)
			#---this is the ion positions after applying nojump
			test = array(ions_traj[ind])
			print 'get distances'
			#---these are the displacements by delta-time (dim 0) and start point (dim 1)
			dists = [[norm(test[i+d,dimslice]-test[i,dimslice])**2 
				for i in range(len(test)-d)] 
				for d in range(len(test))]
			#---the filter is 1 for a valid selection, 0 otherwise, and runs over all frames
			print 'make filter'
			filt = [[0.5*(1*(test[i:i+d+1,2]>zone[0])+1*(test[i:i+d+1,2]<zone[1])) 
				for i in range(len(test)-d)] for d in range(len(test))]
			print 'filter'
			#---we filter the displacements by those with a filter occupancy above a threshold
			dists_filt = [[dists[i][j] for j in range(len(filt[i])) 
				if mean(filt[i][j])>=occupancy] for i in range(len(filt))]
			alldists_filt.append(dists_filt)
		alldists_filt_zone.append(alldists_filt)
	#---get times from the clock
	times = [mean([clock[i+d]-clock[i] for i in range(len(test)-d)]) 
		for d in range(len(test))]

#---plot
if 1:
	fig = plt.figure()
#	ax = plt.subplot(121)
	zdiffusions = [[] for i in range(len(zones))]
	for z in range(len(zones)):
		for l in range(len(alldists_filt_zone[z])):
			# Indexing by zone, then ion, then time point.
			raw = [mean(alldists_filt_zone[z][l][d]) for d in range(500)]
			msd0 =array([[times[k],raw[k]] for k in range(len(times)) if not isnan(raw[k])])
			# msd0 will be empty if the ion is completely absent from the zone for all time
			if len(msd0) > 1:
				ax = plt.subplot(2,round(len(zones)/2.),z+1)
				ax.plot(msd0[:,0],msd0[:,1],c=clrs[z])
				ax.legend()
				fitpart = where((1*(msd0[:,0]>10)+1*(msd0[:,0]<1000))==2)[0]
				[bz,az] = numpy.polyfit(msd0[fitpart,0],msd0[fitpart,1],1)
				zdiffusions[z].append(bz/3/msd_factor)
#		ax.legend((str(round(zones[z][0]-center))+"-"+str(round(zones[z][1]-center))))
		ax.legend(ax,(str(round(zones[z][0]-center))+"-"+str(round(zones[z][1]-center))))
		ax.set_xscale('log')
		ax.set_yscale('log')
	fig = plt.figure()
	ax = plt.subplot(121)
	for z in range(len(zones)):
		ax.plot([z for i in range(len(zdiffusions[z]))],zdiffusions[z],'o',c=clrs[z])
		ax.set_xlim((-1,6))
	ax = plt.subplot(122)
#	ax.plot(range(len(zones)),[mean(zdiffusions[z]) for z in range(len(zones))], 'o-', c='k')
	ax.boxplot([zdiffusions[z] for z in range(len(zones))])
	ax.set_xlim((-1,6))
	# Run show in a separate thread:
	pylab.ion() 
	plt.show()
	

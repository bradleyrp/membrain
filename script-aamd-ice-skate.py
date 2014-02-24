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
		'structure_pkl':'pkl.structures.membrane-v511.a2-surfacer.30000-80000-100.pkl'}}
analysis_names = ['v511-30000-80000-100']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'mset' not in globals():
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset_surf = unpickle(pickles+structure_pkl)
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='cgmd')
		checktime()

# Pseudocode.
# Needs: surface pickle + ion trajectory
# 1. For each ion, from start to end by time (triple loop)
# 		Record the distance from start point by interval length
#		Record if in the zone in separate array (1 or 0)
# 2. Filter array (1 or 0) by some % time in the zone
# 		Only include delta t distances if the mean in the 1/0 array is above some level

#---get ion positions and times
if 1:
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
	ionspos = array(ionspos)[:-1] #---Note: why necessary because 501 frames on one and 500 on other?
	#meshplot(mean(mset_surf.surf[150:]+mean(mset_surf.surf_position),axis=0),vecs=mset_surf.vecs[0])
	#meshpoints(array(ions_traj)[:,0],scale_factor=[20. for i in range(shape(ions_traj)[0])])
	#---NOTE! trajectory[fr].time and trajectory.time don't match! In v510 before I noticed Cal error
#---unwrap PBCs
if 1:
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
	
#---trajectory for a single ion
if 1:
	test = array(ions_traj[0])
	dists = [mean([norm(test[i+d]-test[i])**2 for i in range(len(test)-d)]) for d in range(len(test))]
	times = [mean([clock[i+d]-clock[i] for i in range(len(test)-d)]) for d in range(len(test))]
	plt.plot(times,dists);plt.show()
#---prove pbc
if 1:
	plt.plot(ionspos[:,147,2]);plt.plot(ions_traj[147,:,2]);plt.show()

#---deets
if 1:
	mz = mean(mset_surf.surf_position)
	bufr = [25+25,25+25+20]
	zone1 = [mz+bufr[0],mz+bufr[1]]
	d = 10
	occupancy = 0.1
	nions = len(ions_traj)
	nlips = len(ionspos)
	#---find a better test
	scan = array([[mean(array(ions_traj[i])[:,2]<zone1[1]),mean(array(ions_traj[i])[:,2]>zone1[0])] 
		for i in range(nions)])
	#plt.plot(scan[:,1]+scan[:,0])
	#plt.show()


#---trajectory for a single ion, decomposed
if 1:
	occupancy = 0.5
	incurves = []
	outcurves = []
	allcurves = []
	whichframes = range(nlips)
	for ind in whichframes:
		print ind
		test = array(ions_traj)[ind]
		curv = []
		curv1 = []
		curv1z = []
		curv2 = []
		for d in range(len(ions_traj[0])):
			dists = array([norm(test[i+d][:2]-test[i][:2])**2 for i in range(len(test)-d)])
			inside = (1*(test[:,2]<zone1[1])+1*(test[:,2]>zone1[0])==2)
			outside = (1*(test[:,2]<zone1[1])+1*(test[:,2]>zone1[0])==0)
			allin = [sum(inside[i:i+d+1]) for i in range(len(test)-d)]
			allout = [sum(outside[i:i+d+1]) for i in range(len(test)-d)]
			times = mean([clock[i+d]-clock[i] for i in range(len(test)-d)])
			dists_in = mean(dists[where(array(allin)>d*occupancy)])
			if not isnan(dists_in):
				curv1.append([times,dists_in])
			dists_out = mean(dists[where(array(allin)<d*occupancy)])
			if not isnan(dists_out):
				curv2.append([times,dists_out])
			curv.append([times,mean(dists)])
		curv = array(curv)
		curv1 = array(curv1)
		curv2 = array(curv2)
		incurves.append(curv1)
		outcurves.append(curv2)
		allcurves.append(curv)

#---3D-plotting the not-unwrapped ions
if 0:
	meshplot(mean(mset_surf.surf[150:]+mean(mset_surf.surf_position),axis=0),vecs=mset_surf.vecs[0])
	meshpoints(array(ionspos)[147,:],scale_factor=[10. for i in range(shape(ions_traj)[0])],color=(1,1,1))
	meshpoints(array(ionspos)[414,:],scale_factor=[10. for i in range(shape(ions_traj)[0])],color=(0,1,0))		
	
#---plot
if 1:
	incurves = array(incurves)
	allcurves = array(allcurves)
	inslice = [i for i in range(len(incurves)) if incurves[i] != []]
	allslice = slice(None,None)
	for curv1 in incurves[inslice]:
		if curv1 != []:
			plt.plot(curv1[:,0],curv1[:,1],'r',lw=2,alpha=0.2)
	for curv in allcurves[allslice]:
		plt.plot(curv[:,0],curv[:,1],'b',lw=2,alpha=0.2)
	plt.show()
	

	


#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
analysis_descriptors = {
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl',
		'ionname':'Cal'}}
analysis_names = ['v511-30000-80000-100']
routine = ['compute','postproc'][1:]

#---method
zonesabs = [[0,10],[10,20],[20,30],[40,50],[60,70]]
occupancy = 0.80

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---white whale attempt
if 'compute' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset_surf = unpickle(pickles+structure_pkl)
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---check for pre-existing pickle
		resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		#mdionskate = unpickle(pickles+resultpkl)
		#if mdionskate != None:
		#	raise Exception('except: pkl already exists so figure out your naming problems')
		#---get ion positions and times
		print 'getting ions'
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		vecs=mean(mset_surf.vecs,axis=0)
		ionspos = array(ionspos)[:-1]
		ionstraj = []
		for ind in range(shape(ionspos)[1]):
			course = array(ionspos)[:,ind]
			#---three-line handling PBCs
			hoplistp = (course[1:]-course[:-1])>array(mset_surf.vecs)[1:]/2.
			hoplistn = (course[:-1]-course[1:])>array(mset_surf.vecs)[1:]/2.
			course_nojump = course[1:]-(cumsum(1*hoplistp-1*hoplistn,axis=0))*array(mset_surf.vecs)[1:]
			ionstraj.append(course_nojump)
		ionstraj=array(ionstraj)
		nions = len(ionstraj)
		nframes = len(ionstraj[0])
		center = mean(mset_surf.surf_position)
		thick = mean(mset_surf.surf_thick)
		#---specify zones
		#zones = zones + array(center) + thick
		#zones = list(-1*array(zones)+array(center)-thick)+list(zones+array(center)+thick)
		zones_up = list(zonesabs + array(center) + thick)
		zones_down = list((array(zonesabs)*-1).T[::-1].T + array(center) - thick)
		zones = zones_up + zones_down
		print 'zones'
		print zones
		print vecs
		#---select subset of ions if desired
		ionsel = slice(0,nions)
		nions = len(range(nions)[ionsel])
		#---pre-compute a master array of all displacements
		print 'status: precomputing displacement array, xy'
		dimslice = slice(0,2)
		distsxy = [[norm(ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice],axis=1)**2 
			for i in range(nframes-d)] 
			for d in range(nframes)]
		#---Note, the apply_along_axis code is brutally slow. Better to upgrade to numpy 1.8 on dirac.
		#distsxy = [[np.apply_along_axis(norm,1,ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice])**2 
		#	for i in range(nframes-d)] 
		#	for d in range(nframes)]
		checktime()
		print 'status: precomputing displacement array, z'
		dimslice = slice(2,3)
		distsz = [[norm(ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice],axis=1)**2 
			for i in range(nframes-d)] 
			for d in range(nframes)]
		#distsz = [[np.apply_along_axis(norm,1,ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice])**2 
		#	for i in range(nframes-d)] 
		#	for d in range(nframes)]
		checktime()
		#---loop over zones
		mastermsd_zones = []
		for z in range(len(zones)):
			print 'status: zone = '+str(z)
			zone = zones[z]
			inzone = array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) for i in range(nframes)])
			inzonesliced = [[mean(inzone[i:i+d+1],axis=0) for i in range(nframes-d)] for d in range(nframes)]
			#---scratchwork
			#mask = [np.ma.masked_greater_equal(inzonesliced[d],occupancy) for d in range(nframes)]
			#tmp3 = numpy.ma.masked_greater_equal(inzonesliced[10],0.1)
			#tmp4 = numpy.ma.mask_rowcols(tmp3,axis=1)
			#---scratchwork
			#tmp3 = numpy.ma.masked_greater_equal(inzonesliced[10],0.5)
			#tmp3r = numpy.ma.mask_rows(tmp3)
			#tmp3c = numpy.ma.mask_cols(tmp3)
			#filt = (1*tmp3r.mask+1*tmp3c.mask)
			bigmask = [numpy.ma.masked_greater_equal(inzonesliced[i],occupancy) for i in range(nframes)]
			#filt = [array((1*numpy.ma.mask_rows(bigmask[i]).mask+1*numpy.ma.mask_cols(bigmask[i]).mask)) 
			#	for i in range(nframes)]
			#---Nb this is a fix for masks returning singleton booleans if everything is false
			#mastermsd = [array((1*numpy.ma.mask_rows(bigmask[i]).mask+1*numpy.ma.mask_cols(bigmask[i]).mask)) for i in range(nframes)]
			mastermsd = [(array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) 
				if shape(bigmask[i].mask) != () else zeros(bigmask[i].shape)) for i in range(nframes)]
			mastermsd_zones.append(mastermsd)
			checktime()
		#---record times
		times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
		#---package the results	
		result_data = MembraneData('ionskate')
		#---Nb data are packaged as type (distsxy,distsz,mastermsdzones)
		#---Nb the mastermsd_zones type goes by zone, ion, delta-t, start frame
		result_data.data = [mastermsd_zones,distsxy,distsz]
		result_data.addnote(['times',times])
		result_data.addnote(['occupancy',occupancy])
		result_data.addnote(['zones',zones])
		result_data.addnote(['ionsel',ionsel])
		#---pickle the results
		#pickledump(result_data,resultpkl,directory=pickles)

makeplots = ['panel','calc_diffusion','diffusion'][-1:]

#---postprocess and plot
if 'postproc' in routine:
	#---NOTICE: disabled this to skip writing step because files too big
	if 'skates' not in globals() and 0:
		skates = []
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
			resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
			skate = unpickle(pickles+resultpkl)
			skates.append(skate)
		nframes = 500
		times = array(skate.getnote('times'))
		distsxy = skates[0].get(['type','distsxy'])
		distsz = skates[0].get(['type','distsz'])
	#---select zones for plotting
	#zones_up = list(zonesabs + array(center) + thick)
	#zones_down = list((array(zonesabs)*-1).T[::-1].T + array(center) - thick)
	#zones = zones_up + zones_down
	zones = (zonesabs+[list(i) for i in -1*array(zonesabs)])[:]
	zones_ind = range(len(zones))
	#---you have to tune "upto" according to the time you want, and it's modulated by occupancy
	upto = 6
	times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)][:upto]
	if 'panel' in makeplots:
		print 'status: making panel plots'
		fig = plt.figure()
		gs = gridspec.GridSpec(2,len(zones_ind))
		for z in zones_ind:
			print 'zone = '+str(z)
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[0,z])
			allcurves = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsxy[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clrs[z%len(clrs)],
					markeredgecolor=clrs[z%len(clrs)],alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
		for z in zones_ind:
			print 'zone = '+str(z)
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[1,z])
			allcurves = array([mean(ma.masked_values((1*(z0[i]==2)*distsz[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clrs[z%len(clrs)],alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
		plt.savefig('/home/rpb/ice-panels.png',dpi=300)
		plt.clf()
		plt.close()
	#---calculate diffusion rates
	if 'calc_diffusion' in makeplots:
		print 'status: calculating diffusion rates'
		#---fit to get diffusion rates
		diffusexy = []
		for z in zones_ind:
			print z
			z0 = mastermsd_zones[z]
			allcurves = array([mean(ma.masked_values((z0[i]*distsxy[i]).T,0.).data,axis=1) for i in range(nframes)]).T
			dconsts = []
			for curv in allcurves[::1]:
				curvfilt = array([times,curv[:upto]]).T[curv[:upto]!=0.]
				fitlims = curvfilt[:,0]<1*10**10
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				dconsts.append(bz/3/2)
			diffusexy.append(dconsts)
		diffusez = []
		for z in zones_ind:
			print z
			z0 = mastermsd_zones[z]
			allcurves = array([mean(ma.masked_values((z0[i]*distsz[i]).T,0.).data,axis=1) 
				for i in range(nframes)]).T
			dconsts = []
			for curv in allcurves[::1]:
				curvfilt = array([times,curv[:upto]]).T[curv[:upto]!=0.]
				fitlims = curvfilt[:,0]<1*10**3
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				dconsts.append(bz/3/1)
			diffusez.append(dconsts)
	#---sick plot
	if 'diffusion' in makeplots:
		print 'status: plotting diffusion rates'
		#---specify color gradients for zones
		clrst = [brewer2mpl.get_map('Reds','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		clrsb = [brewer2mpl.get_map('Blues','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		clrsd = clrst+clrsb
		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 40
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in zones_ind[::-1]:
			dat = array([diffusexy[z],diffusez[z]]).T
			ax.plot(dat[:,0],dat[:,1],'o',c=clrsd[z%len(clrsd)],markeredgecolor=clrsd[z%len(clrsd)],alpha=1.0,
				label=str(zones[z][0]/10.)+'-'+str(zones[z][1]/10.)+' nm')
			#---x axis histogram
			counts,edges = histogram(log10(dat[:,0]),bins=nbins)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHistx.plot(bins,counts,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=1.0)
			#---y axis histogram
			bins,counts = histogram(log10(dat[:,1]))
			counts,edges = histogram(log10(dat[:,1]),bins=nbins)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHisty.plot(counts,bins,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=1.0)
		edgeprop = 2.
		xlims = (array(diffusexy).min()/edgeprop,array(diffusexy).max()*edgeprop)
		ylims = (array(diffusez).min()/edgeprop,array(diffusez).max()*edgeprop)
		ax.set_ylim(ylims)
		ax.set_xlim(xlims)
		axHisty.set_ylim(ylims)
		axHistx.set_xlim(xlims)
		ax.set_yscale('log')
		ax.set_xscale('log')		
		ax.grid(True)
		axHistx.grid(True)
		axHisty.grid(True)
		axHistx.set_xscale('log')
		axHisty.set_yscale('log')		
		ax.set_ylabel(r'$D_Z\:(units)$',fontsize=fsaxlabel)
		ax.set_xlabel(r'$D_{XY}\:(units)$',fontsize=fsaxlabel)
		#ax.legend(loc='upper left')
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		axHistx.set_yticklabels([])
		axHisty.set_xticklabels([])
		plt.setp(axHistx.get_xticklabels()+axHisty.get_yticklabels(),visible=False)
		plt.savefig('/home/rpb/ice-diffusions.png',dpi=300)
		#plt.show()
		plt.clf()
		plt.close()
		
#---plotting in 3D with mayavi
if 0:
	meshpoints(ionstraj[72],scale_factor=[5. for i in range(nframes)])
	meshplot(mean(mset_surf.surf,axis=0)+mset_surf.surf_position[0],vecs=mset_surf.vecs[0])

































#---DEVELOPMENT
#-------------------------------------------------------------------------------------------------------------

#---method which loops over zones, but requires too much memory
if 'compute_v2' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset_surf = unpickle(pickles+structure_pkl)
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---check for pre-existing pickle
		resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		mdionskate = unpickle(pickles+resultpkl)
		if mdionskate != None:
			raise Exception('except: pkl already exists so figure out your naming problems')
		#---get ion positions and times
		print 'getting ions'
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		vecs=mean(mset_surf.vecs,axis=0)
		ionspos = array(ionspos)[:-1]
		ionstraj = []
		for ind in range(shape(ionspos)[1]):
			course = array(ionspos)[:,ind]
			#---three-line handling PBCs
			hoplistp = (course[1:]-course[:-1])>array(mset_surf.vecs)[1:]/2.
			hoplistn = (course[:-1]-course[1:])>array(mset_surf.vecs)[1:]/2.
			course_nojump = course[1:]-(cumsum(1*hoplistp-1*hoplistn,axis=0))*\
				array(mset_surf.vecs)[1:]
			ionstraj.append(course_nojump)
		ionstraj=array(ionstraj)
		nions = len(ionstraj)
		nframes = len(ionstraj[0])
		center = mean(mset_surf.surf_position)
		thick = mean(mset_surf.surf_thick)
		#---specify zones
		zones = zones + array(center) + thick
		zones = list(-1*array(zones)+array(center)-thick)+list(zones+array(center)+thick)
		print 'zones'
		print zones
		#---select subset of ions if desired
		ionsel = range(0,nions)
		#---pre-compute a master array of all displacements
		print 'status: precomputing displacement array, xy'
		dimslice = slice(0,2)
		distsxy = [[[norm(ionstraj[k][i+d,dimslice]-ionstraj[k][i,dimslice])**2 for i in 
			range(len(ionstraj[k])-d)] for d in range(len(ionstraj[k]))] for k in ionsel]
		print 'status: precomputing displacement array, z'
		dimslice = slice(2,3)
		distsz = [[[norm(ionstraj[k][i+d,dimslice]-ionstraj[k][i,dimslice])**2 for i in 
			range(len(ionstraj[k])-d)] for d in range(len(ionstraj[k]))] for k in ionsel]
		#---allocate master msd variable
		mastermsd_zones = []
		#---loop over zones
		for z in range(len(zones)):
			zone = zones[z]
			print 'status: computing array of flags for presence in the zones'
			inzone = array([[2==1*(ionstraj[k][i,2]>zone[0])+1*(ionstraj[k][i,2]<zone[1]) 
				for i in range(len(ionstraj[k]))] 
				for k in range(len(ionsel))])	
			print 'status: taking mean occupancy by slicing the zone presence array'
			inzonesliced = array([[[mean(inzone[k][j:j+i+1]) 
				for i in range(nframes-j)] 
				for j in range(nframes)] 
				for k in range(len(ionsel))])
			#print 'status: filtering by occupancy'
			#inzoneslicedocc = array([[[inzonesliced[z][k][j][i]>occupancy 
			#	for i in range(500-j)] 
			#	for j in range(500)] 
			#	for k in range(len(ionsel))])
			#print 'status: applying filter'
			mastermsd = array([[[(distsxy[k][j][i],distsz[k][j][i])
				for i in range(nframes-j) 
				if inzonesliced[k][j][i]>occupancy] 
				for j in range(nframes)] 
				for k in range(len(ionsel))])
			mastermsd_zones.append(mastermsd)
		#---record times
		times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
		#---package the results	
		result_data = MembraneData('ionskate')
		#---Nb data are packaged as type (distsxy,distsz,mastermsdzones)
		#---Nb the mastermsd_zones type goes by zone, ion, delta-t, start frame
		result_data.data = [mastermsd_zones,distsxy,distsz]
		result_data.addnote(['times',times])
		result_data.addnote(['occupancy',occupancy])
		result_data.addnote(['zones',zones])
		result_data.addnote(['ionsel',ionsel])
		#---pickle the results
		pickledump(result_data,resultpkl,directory=pickles)

#---basic plotting outline
if 'plot' in routine:
	fig = plt.figure()
	for z in range(len(zones)):
		ax = plt.subplot(2,round(len(zones)/2.),z+1)
		ax.set_xscale('log')
		ax.set_yscale('log')
		for k in range(len(ionsel)):
			msd = array([[times[j],mean(mastermsd[z][k][j])] for j in range(500) 
				if not isnan(mean(mastermsd[z][k][j]))])
			if msd != []:
				ax.plot(msd[:,0],msd[:,1],'o',c=clrs[z%len(clrs)])
	plt.show()

#---original method works well enough
if 'original' in routine:
	center = mean(mset_surf.surf_position)
	zones = [[0,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]
	zones = zones + array(center)
	zones = zones[:3]
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
	#---pre-compute a master array of all displacements
	print 'status: precomputing displacement array'
	#dists = [[[norm(test[i+d,dimslice]-test[i,dimslice])**2 
	#	for i in range(len(test)-d)] 
	#	for d in range(len(test))] for test in array(ions_traj)[:10]]
	print 'status: computing filter for each zone'
	starttime = time.time()
	for zone in zones:
		alldists_filt = []
		#---deprecated pre-filtering for zones, but mean z is almost always near bilayer
		if 0:
			ions_in_zone = where((1*(meanz>zone[0])+1*(meanz<zone[1]))==2)[0]
			print ions_in_zone
		else:
			ions_in_zone = range(nions)[:10]
		#---compute zone filter once for all ions
		print 'status: precomputing zone filter for all ions'
		filt = [[0.5*(1*(test[i:i+d+1,2]>zone[0])+1*(test[i:i+d+1,2]<zone[1])) 
			for i in range(len(test)-d)] for d in range(len(test))]
		#---loop over ions in the zone and apply the filter
		for ind in ions_in_zone:
			print 'ion = '+str(ind)
			#---this is the ion positions after applying nojump
			#test = array(ions_traj[ind])
			#---these are the displacements by delta-time (dim 0) and start point (dim 1)
			#dists = [[norm(test[i+d,dimslice]-test[i,dimslice])**2 
			#	for i in range(len(test)-d)] 
			#	for d in range(len(test))]
			#---the filter is 1 for a valid selection, 0 otherwise, and runs over all frames
			print 'make filter'
			filt = [[0.5*(1*(test[i:i+d+1,2]>zone[0])+1*(test[i:i+d+1,2]<zone[1])) 
				for i in range(len(test)-d)] for d in range(len(test))]
			print str(time.time()-starttime)
			print 'apply filter'
			#---we filter the displacements by those with a filter occupancy above a threshold
			dists_filt = [[dists[i][j] for j in range(len(filt[i])) 
				if mean(filt[i][j])>=occupancy] for i in range(len(filt))]
			alldists_filt.append(dists_filt)
			alldists_filt_zone.append(alldists_filt)
			print str(time.time()-starttime)
	#---get times from the clock
	times = [mean([clock[i+d]-clock[i] for i in range(len(test)-d)]) for d in range(len(test))]
	
#---original plotting method
if 'original_plot' in routine:
	fig = plt.figure()
	#ax = plt.subplot(121)
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
		#ax.legend((str(round(zones[z][0]-center))+"-"+str(round(zones[z][1]-center))))
		ax.legend(ax,(str(round(zones[z][0]-center))+"-"+str(round(zones[z][1]-center))))
		ax.set_xscale('log')
		ax.set_yscale('log')
	fig = plt.figure()
	ax = plt.subplot(121)
	for z in range(len(zones)):
		ax.plot([z for i in range(len(zdiffusions[z]))],zdiffusions[z],'o',c=clrs[z])
		ax.set_xlim((-1,6))
	ax = plt.subplot(122)
	#ax.plot(range(len(zones)),[mean(zdiffusions[z]) for z in range(len(zones))], 'o-', c='k')
	ax.boxplot([zdiffusions[z] for z in range(len(zones))])
	ax.set_xlim((-1,6))
	# Run show in a separate thread:
	pylab.ion() 
	plt.show()
	
#---vectorized method, but uses too much memory
if 0:
	starttimeall = time.time()
	center = mean(mset_surf.surf_position)
	zones = [[0,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]
	zones = zones + array(center)
	#zones = zones[:3]
	occupancy = 0.8
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
	ionsel = range(0,522,2)
	#---pre-compute a master array of all displacements
	print 'status: precomputing displacement array'
	starttime = time.time()
	dists = [[[norm(ionstraj[k][i+d,dimslice]-ionstraj[k][i,dimslice])**2 for i in 
		range(len(ionstraj[k])-d)] for d in range(len(ionstraj[k]))] for k in ionsel]
	print str(time.time()-starttime)
	#---get inzone
	print 'status: computing array of flags for presence in the zones'
	starttime = time.time()
	inzone = array([[[2==1*(ionstraj[k][i,2]>zone[0])+1*(ionstraj[k][i,2]<zone[1]) 
		for i in range(len(ionstraj[k]))] for k in range(len(ionsel))] for zone in zones])	
	print str(time.time()-starttime)
	#---taking means
	print 'status: taking means'
	starttime = time.time()
	inzonesliced = array([[[[mean(inzone[z][k][j:j+i+1]) for i in range(500-j)] for j in range(500)] 
		for k in range(len(ionsel))] for z in range(len(zones))])
	print str(time.time()-starttime)
	starttime = time.time()
	print 'status: occupancy'
	inzoneslicedocc = array([[[[inzonesliced[z][k][j][i]>0.8 for i in range(500-j)] for j in range(500)] 
		for k in range(len(ionsel))] for z in range(len(zones))])
	print str(time.time()-starttime)
	starttime = time.time()
	print 'status: mastermsd'
	mastermsd = array([[[[dists[k][j][i] for i in range(500-j) if inzonesliced[z][k][j][i]>0.8] 
		for j in range(500)] for k in range(len(ionsel))] for z in range(len(zones))])
	print str(time.time()-starttime)
	print 'done'+str(time.time()-starttimeall)
	times = [mean([clock[i+d]-clock[i] for i in range(500-d)]) for d in range(500)]
	pickledump(mastermsd,'mastermsd-v511',directory=pickles)
	#---plot
	if 0:
		fig = plt.figure()
		for z in range(len(zones)):
			ax = plt.subplot(2,round(len(zones)/2.),z+1)
			ax.set_xscale('log')
			ax.set_yscale('log')
			for k in range(len(ionsel)):
				msd = array([[times[j],mean(mastermsd[z][k][j])] for j in range(500) 
					if not isnan(mean(mastermsd[z][k][j]))])
				if msd != []:
					ax.plot(msd[:,0],msd[:,1],'o',c=clrs[z%len(clrs)])
		plt.show()

#---DEVELOPMENT, JUNK
#-------------------------------------------------------------------------------------------------------------

if 0:
	#occupied = inzonesliced>0.8
	#tmp4 = tmp3[0]*dists
	if 0:
		#---for a single ion and zone
		tmp8 = array([(array(inzonesliced[0][0][i])>0.8) for i in range(500)])
		tmp9 = [[dists[0][d][i]*tmp8[d][i]/sum(tmp8[d][i]) for i in range(500-d)] for d in range(500)]
	#---added if statement
	if 0:
		tmp8 = [1*array([(array(inzonesliced[0][k][i])>0.8) for i in range(500)]) for k in range(len(ionsel))]
		tmp9 = [[[dists[k][d][i]*tmp8[k][d][i]/sum(tmp8[k][d][i]) for i in range(500-d)] 
			for d in range(500)] for k in range(len(ionsel))]
	
	inzone = array([[[2==1*(test[i,2]>zone[0])+1*(test[i,2]<zone[1]) 
		for i in range(len(test))] for test in array(ions_traj)] for zone in zones])
	starttime = time.time()
	inzone = array([[[2==1*(test[i,2]>zone[0])+1*(test[i,2]<zone[1]) 
		for i in range(len(test))] for test in array(ions_traj)] 
		for zone in zones])
	print str(time.time()-starttime)
	starttime = time.time()
	tmp = [[[[mean(inzone[z][k][j:j+i+1]) for i in range(500-j)] 
		for j in range(500)] for k in range(10)] for z in range(len(zones))];
	print str(time.time()-starttime)
	z = 0
	k = 0
	tmp = [[inzone[z][k][j:j+i+1] for i in range(500-j) for j in range(500)]]
	tmp2 = array(tmp)
	for z in range(3):
		tmp4 = tmp3[z]*dists
		plotdat = array([[mean(tmp4[k][i]) for i in range(len(tmp4[k]))] for k in range(10)])
		for pd in plotdat:
			plt.plot(times,pd,c=clrs[z])
	plt.show()
	#filt = [[0.5*(1*(test[i:i+d+1,2]>zone[0])+1*(test[i:i+d+1,2]<zone[1])) 
	#			for i in range(len(test)-d)] for d in range(len(test))]	

	#dists_filt = [[dists[i][j] for j in range(len(filt[i])) 
	#	if mean(filt[i][j])>=occupancy] for i in range(len(filt))]
	
	#---k is ion, i is dt, j is start point
	#dists_filt = [[
	#	dists[k][i][j] 
	#	for j in range(len(dists[i])) 
	#	if mean(inzone[k][j:j+i+1])>=occupancy] 
	#	for i in range(len(dists[k]))] 
	#	for k in range(len(dists))]

	starttime = time.time()
	dists = [[[norm(test[i+d,dimslice]-test[i,dimslice])**2 
		for i in range(len(test)-d)] 
		for d in range(len(test))] for test in array(ions_traj)[:10]]
	print 'done'
	print str(time.time()-starttime)
	

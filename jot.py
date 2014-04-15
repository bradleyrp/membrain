if 0:

	residence_exact = True

	fig=None
	
	normed=True
	scale_by_time=None
	
	zonetype = 'radial_binary'
	
	'''Plot residence time distributions, radial binary version.'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	#---header is the same as the original residence time calculator
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0]+1,[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	jumphists = [histogram(jumpdurations[i],range=(0,ntimes),bins=ntimes+1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	residences = [[] for i in range(len(zonelist))]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[int(jumpbins[k][i])].append(jumpdurations[k][i])
	times = mset.time_dt*histogram(residences[0],range=(1,ntimes+1),bins=ntimes)[1][:-1]
	resdists = [histogram(residences[i],range=(1,ntimes+1),bins=ntimes)[0] 
		for i in range(len(residences)) if i != []]
	zonetype = 'radial_binary'
	times = mset.time_dt*arange(mset.nframes)+10.
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	disctrajt = array(disctraj).T
	t0s = range(0,len(disctrajt)-1)
	t0s = range(0,len(disctrajt)-1)[:1000]
	starttime = time.time()
	init_occupants = []
	for t0 in t0s:
		init_occupants_by_bin = []
		for binn in zonelist:
			init_occupants_by_bin.append(where(disctrajt[t0]==binn)[0])
		init_occupants.append(init_occupants_by_bin)
	collected = [[] for i in zonelist]
	collected_residence = [[] for i in zonelist]
	for t0 in t0s:
		for z in zonelist:
			inds = init_occupants[t0][z]
			staytimes = array([jumptimes[i][jumptimes[i]>t0][0]-t0 for i in inds])
			counts,edges = histogram(staytimes,range=(1,ntimes+1),bins=ntimes)
			collected[z].append(counts)
			if residence_exact:
				residence_curve = [sum(staytimes>i) for i in range(ntimes-t0)]
				collected_residence[z].append(residence_curve)
		status('t0 = '+str(t0),start=starttime,i=t0,looplen=len(t0s))
if 1:
	if residence_exact:
		meanrescurves = [[mean([collected_residence[k][t0][j] 
			for t0 in t0s if j < len(collected_residence[k][t0])]) 
			for j in range(ntimes)] 
			for k in range(len(collected_residence))]
	ax = plt.subplot(111)
	for z in zonelist:
		if normed and not residence_exact:
			#---could normalize by len(mean(collected[z],axis=0)) or sum(mean(collected[z],axis=0)) also
			ax.plot(times,mean(collected[z],axis=0)/mean(collected[z],axis=0)[0],lw=2,
				c=zonecolors[z],label=zonelabels[z])
		elif not normed and not residence_exact:
			ax.plot(times,mean(collected[z],axis=0),lw=2,
				c=zonecolors[z],label=zonelabels[z])
		elif residence_exact:
			ax.plot(times,meanrescurves[z],'.-',lw=1,
				c=zonecolors[z],label=zonelabels[z])
	ax.set_xlim((1,2*times[-1]))
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.grid(True)
	ax.set_ylabel(r'$\mathrm{N_0}$',fontsize=fsaxlabel)
	ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.legend(loc='lower left')
	if residence_exact: extratag = '_dev_exact'
	elif normed: extratag = '_dev_normed'
	else: extratag = '_dev'
	plt.savefig(pickles+'fig-ion_residence'+extratag+'-'+\
		specname_guess(sysname,trajsel).strip('membrane-')+'.png',
		dpi=500,bbox_inches='tight')
	plt.show()


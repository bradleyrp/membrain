
if 0:
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
	starttime = time.time()
	init_occupants = []
	for t0 in t0s:
		init_occupants_by_bin = []
		for binn in zonelist:
			init_occupants_by_bin.append(where(disctrajt[t0]==binn)[0])
		status('t0 = '+str(t0),start=starttime,i=t0,looplen=len(t0s))
		init_occupants.append(init_occupants_by_bin)
	#---new code
	collected = [[] for i in zonelist]
	for t0 in t0s:
		for z in zonelist:
			inds = init_occupants[t0][z]
			staytimes = [jumptimes[i][jumptimes[i]>t0][0]-t0 for i in inds]
			counts,edges = histogram(staytimes,range=(1,ntimes+1),bins=ntimes)
			collected[z].append(counts)
if 1:
	ax = plt.subplot(111)
	for z in zonelist:
		ax.plot(times,mean(collected[z],axis=0),lw=2,
			c=zonecolors[z],label=zonelabels[z])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.grid(True)
	ax.set_ylabel(r'$\mathrm{N_0}$',fontsize=fsaxlabel)
	ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.legend(loc='upper right')
	plt.show()

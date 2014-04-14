



def plot_ion_residence_radial_binary(disctraj,zonetype,fig=None,ignorezero=True,maxtime=100,
	mset=None,dofit=False,cumulative=True):
	'''Plot residence time distributions, .'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0],[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	jumphists = [histogram(jumpdurations[i],range=(0,ntimes),bins=ntimes+1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	#---get all zone definitions for plotting
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	residences = [[] for i in range(len(zonelist))]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[int(jumpbins[k][i])].append(jumpdurations[k][i])
	times = mset.time_dt*histogram(residences[0],range=(1,ntimes),bins=ntimes-1)[1][:-1]
	resdists = [histogram(residences[i],range=(1,ntimes),bins=ntimes-1)[0] 
		for i in range(len(residences)) if i != []]
	#---plot
	ax = plt.subplot(111)
	for r in range(len(resdists)):
		dat = resdists[r]
		c = zonecolors[r]
		alpha = zonealphas[r]
		label = zonelabels[r]
		inds = where(dat!=0)[0]
		ax.plot(times[inds],dat[inds],'o-',c=c,mec=c,lw=1,alpha=alpha,label=label)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((1,ntimes*mset.time_dt*2))
	ax.set_xlabel(r'$\tau\,(ps)$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$frequency$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_title('residence times',fontsize=fsaxtitle)
	ax.legend(loc='upper right')
	ax.grid(True)
	plt.show()


ax = plt.subplot(111)
valid_inds = where(array([len(i) for i in alldat])>0)[0]
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((10**0,max(times)*10))
ax.plot(times[valid_inds],[sum(alldat[i]) for i in valid_inds],'o-')
plt.show()

if 0:
	bin_ind = 3
	ax = plt.subplot(111)
	alldat = [[] for i in range(4999)]
	for ion_ind in range(len(jumpdurations)):
		counts = histogram(jumpdurations[ion_ind][jumpbins[ion_ind]==bin_ind],range=(1,ntimes),bins=ntimes-1)[0]
		for k in where(counts!=0)[0]:
			alldat[k].append(counts[k])
		ax.plot(times[counts!=0],times[counts!=0]*counts[counts!=0]/sum(times[counts!=0]*counts[counts!=0]),'-')
	ax.plot(times,times*mean(alldat,axis=0)/sum(times*mean(alldat,axis=0)),'o-',lw=2,c='k')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((10**0,max(times)*10))
	ax.grid(True)
	ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$\mathrm{R(\mathbf{\tau})\:(ps)}$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_title('time correlations',fontsize=fsaxtitle)

	plt.show()

if 0:
	fig=None
	ignorezero=True
	maxtime=100;
	mset = mset_surf
	dofit=False
	cumulative=False
	'''Plot residence time distributions.'''
	
	spec = 3
	
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0],[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	jumphists = [histogram(jumpdurations[i],range=(1,ntimes),bins=ntimes-1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	residences = [[] for i in range(len(binedges[0])-1)]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[jumpbins[k][i]].append(jumpdurations[k][i])
	times = mset.time_dt*histogram(residences[0],range=(1,ntimes),bins=ntimes-1)[1][:-1]
	resdists = [histogram(residences[i],range=(1,ntimes),bins=ntimes-1)[0] 
		for i in range(len(residences)) if i != []]	

	#---plot
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	halfvec = mean(vecs,axis=0)[2]/2
	mean(vecs,axis=0)[2]
	meanedges = mean(binedges,axis=0)
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	ax = plt.subplot(111)
	for r in range(len(resdists)):
		if ignorezero == False or mids[r] != 0:
			dat = resdists[r]
			if mids[r] > 0.: c = 'r'
			elif mids[r] == 0: c = 'k'
			else: c = 'b'
			ax.plot(times,dat,'o-',c=c,mec=c,lw=1,alpha=1-abs(mids[r]/halfvec),
				label=str(mids[r]/10.)+'('+str(round(-tau,2))+')')
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
	

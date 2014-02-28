#!/usr/bin/python

if 0:
	zone = zones[z]
	print 1
	inzone = array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) for i in range(nframes)])
	print 2
	inzonesliced = [[mean(inzone[i:i+d+1],axis=0) for i in range(nframes-d)] for d in range(nframes)]
	print 3

if 0:
	starttime = time.time()
	bigmask = [numpy.ma.masked_greater_equal(inzonesliced[i],occupancy) for i in range(nframes)]
	print time.time()-starttime
	print 4
	starttime = time.time()
	mastermsd2 = [(array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) if shape(bigmask[40].mask) == () else zeros(bigmask[i].shape)) for i in range(nframes)]
	print time.time()-starttime
	starttime = time.time()
	mastermsd = [array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) for i in range(nframes)]
	print time.time()-starttime

if 0:
	zones_up = list(zonesabs + array(center) + thick)
	zones_down = list((array(zonesabs)*-1).T[::-1].T + array(center) - thick)
	zones = zones_up + zones_down
	inzoneall = [array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) for i in range(nframes)]) for zone in zones]
	#---all MSDs
	times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
	preferred_zones = array([mean(inzoneall[i],axis=0) for i in range(len(zones))]).argmax(axis=0)
if 0:
	clrs = [brewer2mpl.get_map('paired','qualitative',10).mpl_colors[i] for i in range(10)]
	fig = plt.figure()
	ax = plt.subplot(121)
	allrawcurves = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
	for c in range(len(allrawcurves)):
		curv = allrawcurves[c]
		ax.plot(times,curv,c=clrs[preferred_zones[c]])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim((10**-1,10**6))
	ax = plt.subplot(122)
	allrawcurves = array([mean(distsz[i],axis=0) for i in range(nframes)]).T
	for c in range(len(allrawcurves)):
		curv = allrawcurves[c]
		ax.plot(times,curv,c=clrs[preferred_zones[c]])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim((10**-1,10**6))
	plt.savefig('/home/rpb/ice-allmsds.png',dpi=300)
	plt.clf()
	plt.close()

if 1:
	if 0:
		print 'status: calculating diffusion rates'
		#---fit to get diffusion rates
		diffusexy = []
		for z in zones_ind:
			print z
			#print z
			#z0 = mastermsd_zones[z]
			#allcurves = array([mean(ma.masked_values((z0[i]*distsxy[i]).T,0.).data,axis=1) for i in range(nframes)]).T
			allcurves = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T[preferred_zones==z]
			dconsts = []
			for curv in allcurves[::1]:
				curvfilt = array([times,curv]).T
				fitlims = curvfilt[:,0]<1*10**3
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				dconsts.append(bz/3/2)
			diffusexy.append(dconsts)
		diffusez = []
		for z in zones_ind:
			print z
			#print z
			#z0 = mastermsd_zones[z]
			#allcurves = array([mean(ma.masked_values((z0[i]*distsz[i]).T,0.).data,axis=1) for i in range(nframes)]).T
			allcurves = array([mean(distsz[i],axis=0) for i in range(nframes)]).T[preferred_zones==z]
			dconsts = []
			for curv in allcurves[::1]:
				curvfilt = array([times,curv]).T
				fitlims = curvfilt[:,0]<1*10**3
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				dconsts.append(bz/3/1)
			diffusez.append(dconsts)
	if 1:
		print 'status: plotting diffusion rates'
		#---specify color gradients for zones
		clrst = [brewer2mpl.get_map('Reds','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		clrsb = [brewer2mpl.get_map('Blues','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		clrsd = clrst+clrsb
		clrsd = clrs
		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 40
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in zones_ind[::-1]:
			dat = array([diffusexy[z],diffusez[z]]).T
			ax.plot(dat[:,0],dat[:,1],'o',c=clrsd[z%len(clrsd)],markeredgecolor=clrsd[z%len(clrsd)],alpha=0.5,
				label=str(zones[z][0]/10.)+'-'+str(zones[z][1]/10.)+' nm')
			#---x axis histogram
			counts,edges = histogram(log10(dat[:,0]),bins=nbins,normed=True)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHistx.plot(bins,counts,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=0.5)
			#---y axis histogram
			counts,edges = histogram(log10(dat[:,1]),bins=nbins,normed=True)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHisty.plot(counts,bins,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=0.5)
		edgeprop = 2.
		#xlims = (array(diffusexy).min()/edgeprop,array(diffusexy).max()*edgeprop)
		#ylims = (array(diffusez).min()/edgeprop,array(diffusez).max()*edgeprop)
		xlims = (min(min(diffusexy))/edgeprop,max(max(diffusexy))*edgeprop)
		ylims = (min(min(diffusez))/edgeprop,max(max(diffusez))*edgeprop)
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
		plt.savefig('/home/rpb/ice-diffusions-method2.png',dpi=300)
		#plt.show()
		plt.clf()
		plt.close()

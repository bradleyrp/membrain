#!/usr/bin/python

########## original code don't know wtf it is
if 0:
	zone = zones[z]
	print 1
	inzone = array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) 
		for i in range(nframes)])
	print 2
	inzonesliced = [[mean(inzone[i:i+d+1],axis=0) for i in range(nframes-d)] for d in range(nframes)]
	print 3
if 0:
	starttime = time.time()
	bigmask = [numpy.ma.masked_greater_equal(inzonesliced[i],occupancy) for i in range(nframes)]
	print time.time()-starttime
	print 4
	starttime = time.time()
	mastermsd2 = [(array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) 
		if shape(bigmask[40].mask) == () else zeros(bigmask[i].shape)) for i in range(nframes)]
	print time.time()-starttime
	starttime = time.time()
	mastermsd = [array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) for i in range(nframes)]
	print time.time()-starttime

########### preferred zones
if 0:
	zones_up = list(zonesabs + array(center) + thick)
	zones_down = list((array(zonesabs)*-1).T[::-1].T + array(center) - thick)
	zones = zones_up + zones_down
	inzoneall = [array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) 
		for i in range(nframes)]) for zone in zones]
	#---all MSDs
	times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
	preferred_zones = array([mean(inzoneall[i],axis=0) for i in range(len(zones))]).argmax(axis=0)

########## preferred zones, REDUX BITCHES
if 1:
	bwid = 10
	zonesabs = list(arange(center,0-bwid,-bwid))[::-1][:-1]+list(arange(center,vecs[2]+bwid,bwid))
	zonecenters = array(zonesabs[:-1]+zonesabs[1:])/2.
	rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+(array(ionstraj[i,:,2])<0.)*vecs[2] 
		for i in range(nions)]
	rounded = [[int(round(i)) for i in (rewrapz[j]-center)/bwid] for j in range(nions)]
	preferred_zones = stats.mode(rounded,axis=1)[0][:,0]
	rang = int(ptp(preferred_zones))
	pzs = [int(i-preferred_zones.min()) for i in preferred_zones]
	
########### plot allmsds
if 0:
	clrs = [brewer2mpl.get_map('paired','qualitative',12).mpl_colors[i] for i in range(12)]
	fig = plt.figure()
	ax = plt.subplot(121)
	allrawcurves = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
	for c in range(len(allrawcurves)):
		curv = allrawcurves[c]
		ax.plot(times,curv,c=clrs[pzs[c]%len(clrs)])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim((10**-1,10**6))
	ax = plt.subplot(122)
	allrawcurves = array([mean(distsz[i],axis=0) for i in range(nframes)]).T
	for c in range(len(allrawcurves)):
		curv = allrawcurves[c]
		ax.plot(times,curv,c=clrs[pzs[c]%len(clrs)])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim((10**-1,10**6))
	plt.savefig('/home/rpb/ice-allmsds.png',dpi=300)
	plt.clf()
	plt.close()

########### diffusion summary plot of all msds
if 1:
	if 0:
		print 'status: calculating diffusion rates'
		#---fit to get diffusion rates
		allcurvesxy = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
		allcurvesz = array([mean(distsz[i],axis=0) for i in range(nframes)]).T
		diffusexy = []
		for curv in allcurvesxy[::1]:
			curvfilt = array([times,curv]).T
			fitlims = curvfilt[:,0]<1*10**3
			[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
			diffusexy.append(bz/3/2)
		diffusez = []
		for curv in allcurvesz[::1]:
			curvfilt = array([times,curv]).T
			fitlims = curvfilt[:,0]<1*10**3
			[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
			diffusez.append(bz/3/1)
		diffusexy = array(diffusexy)
		diffusez = array(diffusez)
	if 1:
		print 'status: plotting diffusion rates'
		#---specify color gradients for zones
		#clrst = [brewer2mpl.get_map('Reds','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		#clrsb = [brewer2mpl.get_map('Blues','sequential',len(zonesabs)).mpl_colors[i] for i in range(len(zonesabs))][::-1]
		#clrsd = clrst+clrsb
		clrsd = clrs
		#---needs this 
		#allcurvesxy = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
		#allcurvesz = array([mean(distsz[i],axis=0) for i in range(nframes)]).T

		cslice = int(round(center/bwid))
		topinds = range(cslice+1,cslice+rang/2+1)
		botinds = [i%rang for i in range(cslice-1,cslice-1-rang/2,-1)]
		if len(topinds) <= 9:
			redclrs = [brewer2mpl.get_map('Reds','sequential',len(topinds)).mpl_colors[i] for i in range(len(topinds))][::-1]
			bluclrs = [brewer2mpl.get_map('Blues','sequential',len(botinds)).mpl_colors[i] for i in range(len(botinds))][::-1]
		else:
			redclrs = [mpl.cm.RdBu_r(float(i)/rang) for i in range(len(topinds))]
			bluclrs = [mpl.cm.RdBu_r(1.0-float(i)/rang) for i in range(len(botinds))]

		clrsd = [(bluclrs+['k']+redclrs)[(botinds+topinds)[i]] for i in range(len(botinds+topinds))]

		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 40
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in range(rang):
			dat = array([diffusexy[array(pzs)==z],diffusez[array(pzs)==z]]).T	
			if shape(dat) == (2,): dat = array([dat])
			ax.plot(dat[:,0],dat[:,1],'o',c=clrsd[z%len(clrsd)],markeredgecolor=clrsd[z%len(clrsd)],alpha=1.,
				label=str(zonecenters[z]/10.)+' nm')
			#---x axis histogram
			counts,edges = histogram(log10(dat[:,0]),bins=nbins,normed=False)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHistx.plot(bins,counts,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=1.)
			#---y axis histogram
			counts,edges = histogram(log10(dat[:,1]),bins=nbins,normed=False)
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHisty.plot(counts,bins,'o-',lw=2.,c=clrsd[z%len(clrsd)],
				markeredgecolor=clrsd[z%len(clrsd)],alpha=1.)
		edgeprop = 2.
		#xlims = (array(diffusexy).min()/edgeprop,array(diffusexy).max()*edgeprop)
		#ylims = (array(diffusez).min()/edgeprop,array(diffusez).max()*edgeprop)
		xlims = (min(diffusexy)/edgeprop,max(diffusexy)*edgeprop)
		ylims = (min(diffusez)/edgeprop,max(diffusez)*edgeprop)
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
		plt.cla()
		plt.clf()
		plt.close()

if 0:
	mobile_range = (0.02,0.2)
	fig = plt.figure()
	ax = plt.subplot(111)
	for z in range(rang):
		dat = array([diffusexy[array(pzs)==z],diffusez[array(pzs)==z]]).T
		if shape(dat) == (2,): dat = array([dat])
		#---x axis histogram
		counts,edges = histogram(log10(dat[:,0]),bins=10,normed=False,range=[log10(i) for i in mobile_range])
		bins = 10**((edges[:-1]+edges[1:])/2.)
		ax.plot(bins,counts,'o-',lw=2.,c=clrsd[z%len(clrsd)],
			markeredgecolor=clrsd[z%len(clrsd)],alpha=1.)
		ax.set_xscale('log')
		#ax.set_xlim(mobile_range)
	plt.show()

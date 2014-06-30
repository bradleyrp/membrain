#!/usr/bin/python

if 0:
	#---Nb here we analyze points on the surface according to distance to protein (r) and curvature (C0)
	#---method: histogram > smooth > sum frames > norm (r)
	smooth_mean = 1
	#---method: join all C0 for one span, all systems > histogram
	plot_span = 2
	
	#---master list of possible plots
	plot_menu = [
		'hist_sum_norm',
		'hist_smooth_sum_norm',
		'join_systems_histogram',
		'systems_join_span_histogram'
		][-1:]
	for plot_type in plot_menu:
		#---method: histogram > sum frames > norm (r)	
		if plot_type == 'hist_sum_norm':
			fig = plt.figure()
			ax = plt.subplot(111)
			#---sum the histogram over protein distance and C0 magnitude then normalize by protein distance
			dat = sum(bigdathist,axis=0).T/sum(sum(bigdathist,axis=0),axis=1)
			ax.imshow(dat,interpolation='nearest',origin='lower')
			plt.show()
		#---method: histogram > smooth > sum frames > norm (r)	
		if plot_type == 'hist_smooth_sum_norm':
			#---compute
			smoothdat = [scipy.ndimage.filters.gaussian_filter(i,smooth_mean) for i in bigdathist]
			summeddat = sum(smoothdat,axis=0).T/sum(sum(smoothdat,axis=0),axis=1)
			#---plot
			fig = plt.figure()
			ax = plt.subplot(111)
			ax.imshow(summeddat,interpolation='nearest',origin='lower',
				cmap=mpl.cm.jet)
			plt.plot([sum(summeddat[:,i]*range(60))/sum(summeddat[:,i]) for i in range(70)],c='k',lw=3)
			plt.show()
		#---method: for each system join c0 for all frames for all spans > histogram
		if plot_type == 'join_systems_histogram':
			fig = plt.figure()
			ax = plt.subplot(111)
			for j in range(len(md_maps)):
				md_map = md_maps[j]
				for i in range(len(md_map)):
					counts,edges = histogram(array(md_map[i].data).flatten())
					mids = (edges[1:]+edges[:-1])/2.
					plt.plot(mids,counts,lw=2,c=clrs[j])
			plt.show()
		#---method: for each system join c0 for all frames for one span > histogram
		if plot_type == 'systems_join_span_histogram':
			width = 0.35
			nbins = 41
			fig = plt.figure()
			ax = plt.subplot(111)
			ax.axvline(x=0,ymin=0,ymax=1,color='k',lw=2)
			pozsums = []
			negsums = []
			for j in range(len(md_maps)):
				md_map = md_maps[j]
				dat = mean(md_map[plot_span].data,axis=0).flatten()
				lims = [[-1*i,i] for i in [max([abs(dat.min()),abs(dat.max())])]][0]
				counts,edges = histogram(dat,range=lims,bins=nbins)
				mids = (edges[1:]+edges[:-1])/2.
				ax.plot(mids,counts,'o-',lw=2,c=clrs[j],
					label=(analysis_descriptors[analysis_names[j]])['label'],mec=clrs[j])
				pozsums.append(sum(counts[mids>=0]*mids[mids>=0]))
				negsums.append(abs(sum(counts[mids<=0]*mids[mids<=0])))
			#---inset comparing positive and negative values
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="25%",height="40%",loc=1)
			ind = np.arange(len(md_maps))
			pozbars = axins.bar(ind+width,pozsums,width,color='r',alpha=0.5)
			negbars = axins.bar(ind,negsums,width, color='b',alpha=0.5)
			axins.set_xticklabels([(analysis_descriptors[analysis_names[j]])['label'] 
				for j in range(len(analysis_names))])
			plt.setp(axins.get_xticklabels(), rotation=90)
			axins.set_yticklabels([])
			axins.set_xticks(ind+width)
			axins.set_ylabel(r'${C}_{0}\mathbf{(+/-)}$',fontsize=fsaxlabel)
			#---plot
			ax.set_xlabel(r'$\mathsf{C_{0}(x,y)\,(nm^{-1})}$',fontsize=fsaxlabel)
			plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
			ax.set_yticklabels([])
			ax.grid(True)
			ax.legend(loc='upper left')
			plt.savefig(pickles+'fig-stressdist-joinhist-'+bigname+('.flat' if flatz0 else '')+\
				'.png',dpi=300,bbox_inches='tight')
			plt.show()

fig = plt.figure()
ax = plt.subplot(111)
#---sum the histogram over protein distance and C0 magnitude then normalize by protein distance
sourcedat = stack_collected_hists[0]
dat = sum(sourcedat,axis=0).T/sum(sum(sourcedat,axis=0),axis=1)
ax.imshow(dat,interpolation='nearest',origin='lower')
plt.show()


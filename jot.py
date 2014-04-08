if 1:

	example_neighborhood = 6
	peakzoom = [93,108]
	viewplot = False

	if 'descdat' not in globals():
		from scipy import stats
		#---statistics of individual distributions
		drop = 2
		descdat = [[stats.describe(alldats[i][j]) 
			for i in range(shape(alldats)[0])] 
			for j in range(drop,shape(alldats)[1])]
	if 0:
		#---shapiro tests show that most distributions are normal, however some of the large areas are not
		#---this justifies the kurtosis and skew results
		shapiros = array([[scipy.stats.shapiro(x) for x in y] for y in alldats])
	if 1:
		from scipy import stats
		ends = [95,105]
		drop = 2
		statres = array([[stats.ks_2samp(allhists[0][i][ends[0]:ends[1]+1],allhists[2][i][ends[0]:ends[1]+1]),
			stats.ks_2samp(allhists[1][i][ends[0]:ends[1]+1],allhists[2][i][ends[0]:ends[1]+1])]
			for i in range(drop,len(alldats[0]))])
	if 0:
		#---testing different statistics tests
		plt.plot(alldats[0][6]);plt.plot(alldats[1][6],'r');plt.plot(alldats[2][6],'k');plt.show()
		stats.kstest(alldats[0][0],'norm')
		plt.hist(stats.norm.rvs(size=1000));plt.show()
		stats.kstest(stats.t.rvs(3,size=100),'norm')
		plt.hist(stats.t.rvs(2,size=100));plt.show()
		[mean(allhists[i][0]) for i in range(len(allhists))]
		stats.ttest_1samp(allhists[0][0],0.)

	#---master statistics plots
	fig = plt.figure(figsize=(12,12))
	gs1 = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
	gs1.update(top=1.0,bottom=0.7)
	gs2 = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0)
	gs2.update(top=0.62,bottom=0.40)
	gs3 = gridspec.GridSpec(1,4)
	gs3.update(top=0.32,bottom=0.0,wspace=0.25,hspace=0.0)
	
	#---means
	statsaxs = []
	ax = fig.add_subplot(gs3[0])
	tmp = array([[i[2] for i in j] for j in descdat]).T
	for t in range(len(tmp)):
		ax.plot(array(boxslice[drop:])/10.,tmp[t],'o-',c=clrs[t+1],lw=2)
	ax.set_title(r'$\mathrm{mean}\,\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxtitle)
	ax.set_title(r'$\mathrm{mean}$',fontsize=fsaxtitle)
	ax.set_title(r'$\mathrm{mean}\,(nm^{-1})$',fontsize=fsaxtitle)
	if 0: ax.set_ylabel(r'$(nm^{-1})$',fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel,rotation=90)
	statsaxs.append(ax)

	#---variances
	ax = fig.add_subplot(gs3[1])
	tmp = array([[i[3] for i in j] for j in descdat]).T
	for t in range(len(tmp)):
		ax.plot(array(boxslice[drop:])/10.,tmp[t],'o-',c=clrs[t+1],lw=2)
	ax.set_title(r'$\mathrm{variance}\,\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxtitle)
	ax.set_title(r'$\mathrm{variance}$',fontsize=fsaxtitle)
	ax.set_title(r'$\mathrm{variance}\,(nm^{-2})$',fontsize=fsaxtitle)
	if 0: ax.set_ylabel(r'$(nm^{-2})$',fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel,rotation=90)
	statsaxs.append(ax)

	#---skew
	ax = fig.add_subplot(gs3[2])
	tmp = array([[i[4] for i in j] for j in descdat]).T
	for t in range(len(tmp)):
		ax.plot(array(boxslice[drop:])/10.,tmp[t],'o-',c=clrs[t+1],lw=2)
	ax.set_title(r'$\mathrm{skew}\,\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxtitle)
	ax.set_title(r'$\mathrm{skew}$',fontsize=fsaxtitle)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel,rotation=90)
	statsaxs.append(ax)

	#---kurtosis
	ax = fig.add_subplot(gs3[3])
	tmp = array([[i[5] for i in j] for j in descdat]).T
	for t in range(len(tmp)):
		ax.plot(array(boxslice[drop:])/10.,array(tmp[t])+1,'o-',c=clrs[t+1],lw=2)
	ax.set_title(r'$\mathrm{kurtosis}\,\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxtitle)		
	ax.set_title(r'$\mathrm{kurtosis}$',fontsize=fsaxtitle)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel,rotation=90)
	statsaxs.append(ax)
	
	for ax in statsaxs:
		ax.grid(True)
		ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=4))
		ax.set_xlabel(r'$\left|\mathbf{r}_{max}\right|\,\mathrm{(nm)}$',fontsize=fsaxlabel)
		
	#---peak distributions
	ax = fig.add_subplot(gs1[1])
	for j in range(len(allhists)):
		ax.plot(mids[peakzoom[0]:peakzoom[1]+1],allhists[j][6][peakzoom[0]:peakzoom[1]+1],
			'o-',c=clrs[j+1],lw=2)
	ax.grid(True)
	ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=7))
	ax.set_yticklabels([])
	ax.set_xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel)
	ax.set_title('most frequent curvatures',fontsize=fsaxtitle)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)

	#---peak distributions
	ax = fig.add_subplot(gs1[0])
	for j in range(len(allhists)):
		ax.plot(mids,allhists[j][example_neighborhood],'o-',c=clrs[j+1],lw=2,
		label = (analysis_descriptors[analysis_names[j]])['label'])
	ax.legend(loc='upper left')
	ax.grid(True)
	ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=7))
	ax.set_xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel)
	ax.set_yticklabels([])
	ax.set_title('curvatures',fontsize=fsaxtitle)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)

	#---C0 balances
	lims = [10**10,0]
	axs = []
	for i in range(len(allpoz)):
		ax = fig.add_subplot(gs2[i])
		ax.plot(array(rvals)/10.,allpoz[i][-len(rvals):],'o-',
			label=r'$\left\langle C_{0}>0 \right\rangle$ ',
			color=clrs[i+1])
		ax.plot(array(rvals)/10.,allneg[i][-len(rvals):],'s-',
			label=r'$\left\langle C_{0}<0 \right\rangle$ ',
			color=clrs[i+1])
		ax.set_title((analysis_descriptors[analysis_names[i]])['label'],fontsize=fsaxtitle)
		lims[1] = max([max(allneg[i]),max(allpoz[i])]) \
			if max([max(allneg[i]),max(allpoz[i])]) > lims[1] else lims[1]
		lims[0] = min([min(allneg[i]),min(allpoz[i])]) \
			if min([min(allneg[i]),min(allpoz[i])]) < lims[0] else lims[0]
		if i > 0: ax.set_yticklabels([])
		if i == 0: 
			if 0: ax.set_ylabel(r'$\left\langle C_{0} \right\rangle (\mathrm{{nm}^{-1}})$',fontsize=fsaxlabel)
			ax.set_ylabel(r'$\mathrm{C_0\,balance}$',fontsize=fsaxlabel)
		axs.append(ax)
	for ax in axs:
		ax.legend(loc='lower right')
		ax.set_ylim((lims[0]*0.95,lims[1]*1.01))
		ax.grid(True)
		ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=7))
		ax.set_xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel)
		ax.set_yticklabels([])
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		ax.set_xlabel(r'$\left|\mathbf{r}_{max}\right|\,\mathrm{(nm)}$',fontsize=fsaxlabel)
	
	plt.savefig(pickles+'fig-stressmap-overview-'+bigname+('.flat' if flatz0 else '')+\
		'-span'+str(span_sweep[plot_span])+\
		('-blur'+str(smooth_wid) if do_blur else '')+\
		('-blurfirst'+str(smooth_wid) if do_blur_first else '')+\
		'.png',dpi=300,bbox_inches='tight')

	if viewplot: plt.show()
	plt.close()
	print 'done'



if 0:
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		msdat = unpickle_special('pkl.dimple2avg.'+filespec+'.pkl')
		anum = analysis_names.index(aname)
		msdats[anum] = msdat
	if 0:
		hmaxs_by_cutoff = []
		for param in params_plot:
			cutoff,zfiltdir = param
			#---lookup correct msdat items		
			msdats_subset = [[] for i in range(len(analysis_names))]
			for aname in analysis_names:
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
				anum = analysis_names.index(aname)
				for ms in range(len(msdats[anum])):
					if msdats[anum][ms].getnote('cutoff') == cutoff and \
						msdats[anum][ms].getnote('zfiltdir') == zfiltdir and\
						msdats[anum][ms].getnote('decay_z0_min') == None:
						msdats_subset[anum].append(ms)
			hmaxs_by_cutoff.append([[msdats[j][i].get(['type','maxhs'])[0] 
				for i in msdats_subset[j]] for j in range(len(msdats))])
	cuts = [i[0] for i in params_plot]
	cutoff = cuts[0]
	zfiltdir = 0
	msdats_subset = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		for ms in range(len(msdats[anum])):
			if msdats[anum][ms].getnote('cutoff') == cutoff and \
				msdats[anum][ms].getnote('zfiltdir') == zfiltdir and \
				msdats[anum][ms].getnote('decay_z0_min') == None:
				msdats_subset[anum].append(ms)
	cuts = [i[0] for i in params_plot]
	cutoff = cuts[0]
	zfiltdir = 0
	msdats_subset2 = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		for ms in range(len(msdats[anum])):
			if msdats[anum][ms].getnote('cutoff') == cutoff and \
				msdats[anum][ms].getnote('zfiltdir') == zfiltdir and \
				msdats[anum][ms].getnote('decay_z0_min') == True:
				msdats_subset2[anum].append(ms)
if 0:
	fig = plt.figure(figsize=(10,len(analysis_names)*2))
	gs = gridspec.GridSpec(len(analysis_names),3,wspace=0.0,hspace=0.0)
	cuts = [i[0] for i in params_plot]	
	extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())]) for mset in msets])
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		mset = msets[anum]
		#---get protein
		if protein_pkl == None: mset_protein = msets[anum]
		#---retrieve protein points from a different dictionary item
		else:
			protein_pkl_name = protein_pkl
			for i in analysis_descriptors[protein_pkl_name]: 
				vars()[i] = (analysis_descriptors[protein_pkl_name])[i]
			filespec = specname_guess(sysname,trajsel)
			mset_protein = unpickle_special('pkl.structures.'+filespec+'.pkl')
			#---re-populate the variables for the original system now that you got the protein pickle
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---protein colors
		dat = [msdats[anum][i] for i in msdats_subset[anum]]
		casual_names = [test.getnote('testlist')[test.getnote('this_test')] for test in dat]
		casual_names = [i[0] if type(i) == list else i for i in casual_names]
		color_order = [casual_names.index(i) for i in casual_names if i not in ['peak','valley','full']]
		#---plot Hmax values for the average structures
		ax = plt.subplot(gs[anum,0])
		testnames = [i[0] if type(i) == list else i for i in msdats[anum][0].getnote('testlist')]
		print len(dat)
		print testnames
		print casual_names
		print color_order
		for testnum in range(len(testnames)):
			test = dat[testnum]
			testcurv = [[i.get(['type','maxhs'])[0] for i in msdats[anum]
				if (i.getnote('cutoff') == cutoff and \
				i.getnote('decay_z0_min') == None and \
				i.getnote('this_test') == testnum)][0] 
				for cutoff in cuts]
			casual_name = test.getnote('testlist')[test.getnote('this_test')]
			if casual_name == 'peak': color = colordict('blue')
			elif casual_name == 'valley': color = colordict('red')
			elif casual_name == 'full': color = colordict('black')
			else: color = colordict(color_order.index(testnum))
			ax.plot(array(cuts)/10.,testcurv,'-',c=color,label=label,lw=2)
		ax.axhline(y=0,xmax=1,xmin=0,lw=2,c='k')
		ax.set_ylim((-0.02,0.02))
		ax.grid(True)
		if analysis_names.index(aname) < len(analysis_names)-1: ax.set_xticklabels([])
		else: ax.set_xlabel(r'cutoff (nm)',fontsize=fsaxlabel)
		ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
		if analysis_names.index(aname) == 0:
			ax.set_title('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
			
		#---plot Hmax values for the average structures, z0_min
		ax = plt.subplot(gs[anum,2])
		testnames = [i[0] if type(i) == list else i for i in msdats[anum][0].getnote('testlist')]
		print len(dat)
		print testnames
		print casual_names
		print color_order
		for testnum in range(len(testnames)):
			test = dat[testnum]
			testcurv = [[i.get(['type','maxhs'])[0] for i in msdats[anum] 
				if (i.getnote('cutoff') == cutoff and \
				i.getnote('decay_z0_min') == True and \
				i.getnote('this_test') == testnum)][0] 
				for cutoff in cuts]
			casual_name = test.getnote('testlist')[test.getnote('this_test')]
			if casual_name == 'peak': color = colordict('blue')
			elif casual_name == 'valley': color = colordict('red')
			elif casual_name == 'full': color = colordict('black')
			else: color = colordict(color_order.index(testnum))
			ax.plot(array(cuts)/10.,testcurv,'-',c=color,label=label,lw=2)
		ax.axhline(y=0,xmax=1,xmin=0,lw=2,c='k')
		ax.set_ylim((-0.02,0.02))
		ax.grid(True)
		if analysis_names.index(aname) < len(analysis_names)-1: ax.set_xticklabels([])
		else: ax.set_xlabel(r'cutoff (nm)',fontsize=fsaxlabel)
		ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
		if analysis_names.index(aname) == 0:
			ax.set_title('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
		ax.set_ylabel(label,fontsize=fsaxlabel)
		#---plot structure
		ax = plt.subplot(gs[anum,1])
		im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=mpl.cm.RdBu_r,
			inset=False,cmap_washout=1.0,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.])
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_ylabel('')
		ax.set_xlabel('')
		#---height color scale
		axins2 = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins2,orientation="vertical")
		axins2.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
		if analysis_names.index(aname) == 0:
			ax.set_title(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$')
		#---protein hulls
		for pnum in range(len(dat)):
			test = dat[pnum]
			ps = test.getnote('tl')[test.getnote('this_test')]
			if type(ps) == slice:
				protpts = mean(mset_protein.protein,axis=0)[ps,0:2]
			else: 
				protpts = mean(mset_protein.protein,axis=0)[:,0:2]+ps
			casual_name = test.getnote('testlist')[test.getnote('this_test')]
			if casual_name != 'full' or len(casual_names) == 1:
				if casual_name == 'peak': color = colordict('blue')
				elif casual_name == 'valley': color = colordict('red')
				elif casual_name == 'full': color = colordict('black')
				else: color = colordict(color_order.index(pnum))
				plothull(ax,protpts,mset=mset_protein,subdivide=None,c=color,alpha=1.)
		fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
		plt.savefig(pickles+'fig-dimple2avg-'+bigname+'.png',dpi=300,bbox_inches='tight')
		if show_plots: plt.show()
		plt.clf()
if 0:
	mset = unpickle(pickles+'pkl.structures.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl')
	plotter_undulate(mset)
if 0:
	mset = unpickle(pickles+'pkl.cells.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl')


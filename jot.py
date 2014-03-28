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

mset = unpickle(pickles+'pkl.structures.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl')
plotter_undulate(mset)

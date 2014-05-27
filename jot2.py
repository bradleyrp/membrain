#!/usr/bin/python

def plot_ion_distribution(disctraj_prime,binedges,fig=None,thisax=None,text=True,
	colors=True,mset=None,barplot=False,bintype=None,ionlist=None,ionname=None,ls=None,label=None,
	bulk_relative=True):
	'''Generic plotting function for ion distributions by height from midplane according to binedges.'''
	peakval = 0.
	#---Nb only works with discretizer_z
	if fig == None: fig = plt.figure()
	if mset == None: mset = globals()['mset']
	gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
	if thisax == None: ax = fig.add_subplot(gs[0])
	else: ax = thisax
	counts,edges = histogram(disctraj.flatten(),range=(0,len(binedges[0])-1),bins=len(binedges[0])-1,
		normed=True)
	meanedges = mean(array(binedges),axis=0)
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	lefts = meanedges[:-1]	
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	mids = [round(i) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	halfvec = mean(vecs,axis=0)[2]/2
	if bintype == 'fixed': zslice = slice(1,-1)
	else: zslice = slice(None,None)
	if bulk_relative:
		#---buffer in Angstroms at maximum distance from bilayer to compute bulk
		buffer_size = 20
		meanbins = mean(binedges,axis=0)
		topbin = where(meanbins>(meanbins.max()-buffer_size))[0][0]
		botbin = where(meanbins<(-meanbins.max()+buffer_size))[0][-1]
		bulkconc = mean([mean(counts[1:botbin]),mean(counts[topbin:-1])])
		print 'status: bulk concentration = '+str(bulkconc)
	else: bulkconc = 1.
	if barplot:
		for n in (range(len(mids))[slice(1,len(mids-1))] if bintype == 'fixed' else range(len(mids))):
			if mids[n] == 0: color = 'k'
			elif mids[n] > 0: color = 'r'
			elif mids[n] < 0: color = 'b'
			alpha = 1-abs(mids[n]/max([halfvec,max(mids)]))
			ax.bar(lefts[n],counts[n],color=color,width=list(meanedges[1:]-meanedges[:-1])[n],
				alpha=1)
	else:
		#---color selection section
		if resname_group == 'phosphate_position':
			color = color_dictionary_aamd(ionname=ion_name,lipid_resname=ptdins_resname,
				comparison='ions_phospate_position')
		elif resname_group == 'protonation':
			color = color_dictionary_aamd(ionname=ion_name,
				lipid_resname=ptdins_resname,comparison='protonation')
		else:
			color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
		#---label definitions
		if resname_group == 'phosphate_position':
			label = proper_ion_labels[ionlist[ionnum]]+' with '+\
				(proper_ion_labels[ionlist[0]]+' and ' if ionnum>0 else '')+\
				extra_label_list[anum]
		elif resname_group == 'protonation':
			label = proper_ion_labels[ionlist[ionnum]]+' with '+extra_label_list[anum]
		else:
			label=label
		#---plot
		ax.plot(1./2*(meanedges[1:]+meanedges[:-1])[zslice],counts[zslice]/bulkconc,ls,
			color=color,
			alpha=1,lw=2.5,label=label)
		if peakval < array(counts[zslice]/bulkconc).max(): peakval = array(counts[zslice]/bulkconc).max()
	ax.grid(True)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_xlabel(r'$\mathrm{z(x,y)\,\mathrm{(\AA)}}$',fontsize=fsaxlabel)
	ax.set_ylabel(('relative concentration' if bulk_relative == True else 'concentration'),
		fontsize=fsaxlabel)
	ax.set_xlim((min(1./2*(meanedges[1:]+meanedges[:-1])),max(1./2*(meanedges[1:]+meanedges[:-1]))))
	ax.set_yticks(list(arange(0,1.1*peakval,1))[:-1])
	ax.set_yticklabels(list(arange(0,1.1*peakval,1))[:-1])
	if len(list(arange(0,1.1*peakval,1))[:-1]) > 10:
		ax.set_yticks([1]+list(arange(4,1.1*peakval,4)))
		ax.set_yticklabels([1]+list(arange(4,1.1*peakval,4)))
	ax.set_ylim((0,1.1*peakval))
	if text == False:
		ax.set_ylabel('')
		ax.set_xlabel('')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
	if fig == None: plt.show()
	return array(counts[zslice]/bulkconc),edges

#---example for doing the coordinate shift and the binning in the z-dimension
if 'compute_z' in routine:
	fig = plt.figure(figsize=(10,6))
	gs = gridspec.GridSpec(2,1,wspace=0.0,hspace=0.0)
	ax = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1])
	axes = [ax,ax2]
	lslist = ['-','-']
	peakval = [0,0]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		mset = msets[anum]
		mset_surf = msets_surf[anum]
		ionlist = master_ionlist[anum]
		binedges,monobins = binlister('fixed',mset=mset_surf,binw=1,monoz=master_monoz[anum])
		valid_inds = where([len(i)==int(round(mean([len(i) 
			for i in binedges]))) for i in binedges])[0]
		binedges = list(array(binedges)[valid_inds])
		for ionnum in range(len(ionlist)):
			relpos = bilayer_shift(mset=mset,vecs=array([mset.vec(i) for i in range(mset.nframes)]),
				ionspos=master_ionspos[anum][ionnum],midz=master_midz[anum])
			disctraj = discretizer_z(binedges,array(relpos)[valid_inds])
			counts,edges = plot_ion_distribution(disctraj,binedges,
				ionname=ionlist[0],ionlist=ionlist,
				label=proper_ion_labels[ionlist[ionnum]]+\
					(' with '+proper_ion_labels[ionlist[0]] if ionnum>0 else ''),
				bintype='fixed',fig=fig,thisax=axes[ionnum],ls=lslist[ionnum],
				bulk_relative=norm_z_concentration)
			if peakval[ionnum] < counts.max():
				peakval[ionnum] = counts.max()
	ax.legend(loc='upper left',fontsize=fsaxlegend)
	ax2.legend(loc='upper left',fontsize=fsaxlegend)

	#---note that the following plot specs were originally in the function, but I moved them here
	if peakval[1] < 3: peakval[1] = 3
	ax.set_yticks(list(arange(0,1.1*peakval[0],1))[:-1])
	ax.set_yticklabels([str(int(i)) for i in list(arange(0,1.1*peakval[0],1))[:-1]])
	if len(list(arange(0,1.1*peakval[0],1))[:-1]) > 10:
		ax.set_yticks([1]+list(arange(4,1.1*peakval[0],4)))
		ax.set_yticklabels([1]+list([str(int(i)) for i in  arange(4,1.1*peakval[0],4)]))
	ax.set_ylim((0,1.1*peakval[0]))
	ax2.set_yticks(list(arange(0,1.1*peakval[1],1))[:-1])
	ax2.set_yticklabels([str(int(i)) for i in list(arange(0,1.1*peakval[1],1))[:-1]])
	if len(list(arange(0,1.1*peakval[1],1))[:-1]) > 10:
		ax2.set_yticks([1]+list(arange(4,1.1*peakval[1],4)))
		ax2.set_yticklabels([1]+list([str(int(i)) for i in  arange(4,1.1*peakval[1],4)]))
	ax2.set_ylim((0,1.1*peakval[1]))

	ax.set_xlabel('')
	ax.set_xticklabels([])
	ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
	#---save
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-ion_equilibrium_z-'+\
		('not_normalized-' if norm_z_concentration == False else '')+\
		'-'.join(analysis_names)+'.png',dpi=300)
	if showplots: plt.show()
	plt.close(fig)


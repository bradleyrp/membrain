#!/usr/bin/python

logfile,interact,debugmode = [None,False,None]
from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

if 'batch_override' not in globals():

	#---parameters
	rounder = 4

	#---settings
	compare_phosphate_position = False

	#---standard selection
	analysis_names = [
		'v530-40000-90000-50',
		'v531-40000-90000-50',
		'v532-40000-90000-50',
		'v509-40000-90000-50',
		'v510-40000-90000-50',
		'v511-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		'v514-22000-32000-10',
		'v515-20000-30000-10',
		][:1]
	
	#---alternate tests
	if compare_phosphate_position:
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	
	#---routine
	routine = [
		'load',
		'compute',
		][:1]
	
	#---settings
	showplots = False

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'load' in routine:
	#---master data
	master_sortdat = []
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		status('status: system = '+aname+'\n')
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		traj = trajfile[0]

		#---note: removed a check to see if the pkl exists from this space
		mset = MembraneSet()
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		mset.identify_monolayers(director)
		if residues == 'infer':
			residues = list(set(mset.universe.residues.resnames()))
		mset.identify_residues(residues)
		msets.append(mset)
		
		#---load all of the residues separately by monolayer
		bigsel = mset.universe.selectAtoms(selector)
		nmonos = 2
		nrestypes = len(mset.resnames)
		if whichframes == None: frameslice = range(len(mset.universe.trajectory))
		else: frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
		sortdat = [[[] for i in range(nrestypes)] for j in range(nmonos)]
		for frameno in frameslice:
			status('status: frame = '+str(frameno+1)+'/'+str(len(frameslice)))
			mset.gotoframe(frameno)
			mset.vec(frameno)
			for monon in range(nmonos):
				for resn in range(nrestypes):
					if mset.monolayer_by_resid[monon][resn] != []:
						sortdat[monon][resn].append(bigsel.coordinates()[
							array(mset.monolayer_by_resid[monon][resn])])
					else: sortdat[monon][resn].append([])
		status('status: done'+'\n')
		master_sortdat.append(sortdat)
				

#---LOAD
#-------------------------------------------------------------------------------------------------------------

if 'compute' in routine:

	msds_master = [[] for aname in analysis_names]
	time_master = []
	reslist = [[] for aname in analysis_names]
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		mset = msets[anum]
		dists_all = [[[] for i in range(nrestypes)] for j in range(nmonos)]
		for resn in range(shape(master_sortdat[anum])[1]):
			resname = mset.resnames[resn]
			tmp = []
			for monon in range(2):
				if shape(master_sortdat[anum][monon][mset.resnames.index(resname)])[1] != 0:
					tmp.append(array(master_sortdat[anum][monon]\
						[mset.resnames.index(resname)]).transpose(1,0,2))
			if len(shape(tmp)) == 4 and shape(tmp)[0] == 1: trajec = tmp[0]
			elif len(shape(tmp)) == 4 and shape(tmp)[0] > 1: 
				trajec = concatenate((tmp[0],tmp[1]),axis=0)
			elif len(shape(tmp)) == 1:
				print 'warning: you have a lipid flip with residue counts: '+\
					str(len(tmp[0]))+','+str(len(tmp[1]))+' but I will procede'
				trajec = concatenate((tmp[0],tmp[1]),axis=0)
			else: trajec = []
			if trajec != []:
				reslist.append(resname)
				dimslice = slice(0,3)
				dtlimit = shape(trajec)[1]
				nframes = shape(trajec)[1]
				sel = array(range(shape(trajec)[0]))
				distsxyz = [[norm(trajec[sel,i+d,dimslice]-trajec[sel,i,dimslice],axis=1)**2/100. 
					for i in range(nframes-d)] 
					for d in range(dtlimit)]
				msds = array([mean(distsxyz[i],axis=0) for i in range(dtlimit)]).T	
				times = mset.time_dt*array(range(shape(trajec)[1]))
				msds_master[anum].append(array([mean(distsxyz[i],axis=0) for i in range(dtlimit)]).T)
				status('status: computing MSD, resname = '+resname+'\n')
				del distsxyz
		time_master.append(times)


if 'impressionist' in routine:

	#---colorful plot of lipid motions
	back_lipid = array(sortdat[0][mset.resnames.index('DOPE')])
	front_lipid = array(sortdat[0][mset.resnames.index('PI2P')])
	vecs = mean(mset.vecs,axis=0)
	fig = plt.figure()
	ax = plt.subplot(111)
	#---plot colorful background lipids
	
	for resi in range(len(back_lipid[0])):
		for shift in [[i,j] for i in [-1,0,1] for j in [-1,0,1]]:
			plt.plot((back_lipid[:,resi,0]+shift[0]*vecs[0])/10.,
				(back_lipid[:,resi,1]+shift[1]*vecs[1])/10.,c=np.random.rand(3,1),
				alpha=0.5)
	#---box vectors behind foreground lipid
	size_x = abs(-vecs[0]/2.-vecs[0]*1.5)
	size_y = abs(-vecs[1]/2.-vecs[1]*1.5)
	ax.set_xlim((-vecs[0]/2./10.,vecs[0]*1.5/10.))
	ax.set_ylim((-vecs[1]/2./10.,vecs[1]*1.5/10.))
	ax.axvline(x=0,ymin=0.25,ymax=0.75,color='k',lw=2)
	ax.axvline(x=vecs[0]/10.,ymin=0.25,ymax=0.75,color='k',lw=2)
	ax.axhline(y=0,xmin=0.25,xmax=0.75,color='k',lw=2)
	ax.axhline(y=vecs[1]/10.,xmin=0.25,xmax=0.75,color='k',lw=2)
	#---plot white foreground lipid
	for resi in range(len(front_lipid[0])):
		for shift in [[i,j] for i in [-1,0,1] for j in [-1,0,1]]:
			plt.plot((front_lipid[:,resi,0]+shift[0]*vecs[0])/10.,
				(front_lipid[:,resi,1]+shift[1]*vecs[1])/10.,c='r')
	ax.set_aspect(1)
	ax.set_title(ptdins_label+' motion',fontsize=fsaxlabel)
	ax.set_xlabel(r'$\mathrm{X\:(nm)}$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$\mathrm{Y\:(nm)}$',fontsize=fsaxlabel)
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-impressionist-'+aname+'.png',dpi=300)
	if showplots: plt.show()
	plt.close(fig)

if 'plot_lipid_diffusion_detail' in routine:

	#---note the following section just handles different naming conventions, it works but it's ugly
	all_resnames = list(set(concatenate([i.resnames for i in msets])))
	all_resnames = [('ptdins' if i in key_ptdins else i) for i in all_resnames]
	valid_resnames = [r for r in all_resnames if 
		all([
			((analysis_descriptors[analysis_names[msets.index(mset)]])['ptdins_resname']
			if r == 'ptdins' else r) 
			in mset.resnames for mset in msets])]
	valid_resnames = list(set(valid_resnames))

	for resname in valid_resnames:

		#---figure
		fig = plt.figure(figsize=(8,3*len(analysis_names)))
		gs = gridspec.GridSpec(len(analysis_names),1,wspace=0.0,hspace=0.0)

		#---settings
		show_all = True
		fix_lims = [(100,100000),(10**-2,10**1)]

		#---loop over analysis questions
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			mset = msets[anum]
			if resname == 'ptdins': resname_actual = ptdins_resname
			else: resname_actual = resname
			ax = fig.add_subplot(gs[anum])
			
			#---fit
			msds_mean = mean(msds_master[anum][mset.resnames.index(resname_actual)],axis=0)
			datbuf = [int(len(msds_mean)*.1)+2,int(len(msds_mean)*.9)+1]
			[q,r], covariance = polyfit(time_master[anum][datbuf[0]:datbuf[1]],
				msds_mean[datbuf[0]:datbuf[1]],1,cov=True)
			print 'aname = '+aname.ljust(22)+'resname = '+resname.ljust(6)+' D = ' + '{:0.3f}'.format(q*10**6/4.).ljust(10) +' +/- ' + '{:0.3f}'.format(sqrt(covariance[0][0])*10**6/4).ljust(8) + ' (micron**2/second)'
			
			#---plot all curves
			if show_all:
				#---color selection section
				if resname_group == 'phosphate_position':
					color = color_dictionary_aamd(ionname=ion_name,lipid_resname=ptdins_resname,
						comparison='ions_phospate_position')
				elif resname_group == 'protonation':
					color = color_dictionary_aamd(ionname=ion_name,
						lipid_resname=ptdins_resname,comparison='protonation')
				else:
					color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
				for m in msds_master[anum][mset.resnames.index(resname_actual)]:
					ax.plot(time_master[anum],m,alpha=0.35,lw=2,color=color)
				ax.plot(time_master[anum],msds_mean,'k--',lw=2)
				#---label definitions
				if resname_group == 'phosphate_position':
					label=ion_label+extra_label_list[anum]+'\n'+\
						'D = '+'{0:.2f}'.format(q*10**6/4.)+\
						r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
				elif resname_group == 'protonation':
					if resname_actual != ptdins_resname:
						label = 'with '+ion_label+' and '+\
							extra_label_list[anum]+'\n'+\
								'D = '+'{0:.2f}'.format(q*10**6/4.)+\
								r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
					else:
						label = extra_label_list[anum]+' and '+ion_label+\
							'\n'+'D = '+'{0:.2f}'.format(q*10**6/4.)+\
								r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
				else:
					label=ion_label+', D = '+'{0:.2f}'.format(q*10**6/4.)+\
						r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
				p1 = plot_with_rectangle(time_master[anum][datbuf[0]:datbuf[1]],
					msds_mean[datbuf[0]:datbuf[1]],
					ax=ax,new_color=color,color='k',lw=2,
					label=label)
				handles, labels = ax.get_legend_handles_labels()
				ax.legend(handles[1:],labels[1:],loc='upper left')
			#---plot settings
			ax.set_xscale('log')
			ax.set_yscale('log')
			if fix_lims != None:
				ax.set_xlim(fix_lims[0])
				ax.set_ylim(fix_lims[1])
			if 0: ax.grid(b=True,which='major',color='k',linestyle='-')
			if 0: ax.grid(b=True,which='minor',color='k',linestyle='--',alpha=0.2)
			ax.grid(b=True,which='major',color='k',linestyle='--',alpha=0.5)
			ax.set_xlabel('Time (ps)',fontsize=fsaxlabel)
			ax.set_ylabel('MSD (nm$^{2}$)',fontsize=fsaxlabel)
			plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
			plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper'))
			if anum < len(analysis_names)-1:
				ax.set_xticklabels([])
				ax.set_xlabel('')
			if anum == 0: 
				ax.set_title('Lateral diffusion, '+\
					(proper_residue_names_long[resname_actual] if resname == resname_actual
					else 'PtdIns'),fontsize=fsaxlabel)

		#---save
		plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-lipid-diffusion-'+resname+'-'+\
			'-'.join(analysis_names)+'.png',dpi=300)
		if showplots: plt.show()
		plt.close(fig)

if 'plot_lipid_diffusion_summary' in routine:

		#---note the following section just handles different naming conventions, it works but it's ugly
		all_resnames = list(set(concatenate([i.resnames for i in msets])))
		all_resnames = [('ptdins' if i in key_ptdins else i) for i in all_resnames]
		valid_resnames = [r for r in all_resnames if 
			all([
				((analysis_descriptors[analysis_names[msets.index(mset)]])['ptdins_resname']
				if r == 'ptdins' else r) 
				in mset.resnames for mset in msets])]
		valid_resnames = list(set(valid_resnames))
		#---ptdins goes first
		if any([i in valid_resnames for i in key_ptdins]) or 'ptdins' in valid_resnames:
			if 'ptdins' in valid_resnames:
				valid_resnames = ['ptdins']+[i for i in valid_resnames if i != 'ptdins']
			else: 
				valid_resnames = [ptdins_resname]+[i for i in valid_resnames if i not in key_ptdins]
		#---figure
		fig = plt.figure(figsize=(8,2.5*len(valid_resnames)))
		gs = gridspec.GridSpec(len(valid_resnames),1,wspace=0.0,hspace=0.0)
		#---settings
		show_means = True
		fix_lims = [(5*10**2,100000),(10**-1,2*10**0)]
		#---temporary setting
		if resname_group == 'protonation': fix_lims = [(5*10**1,100000),(10**-1,2*10**0)]
		for resname in valid_resnames:
			pnum = valid_resnames.index(resname)
			ax = fig.add_subplot(gs[pnum])
			#---loop over analysis questions
			for aname in analysis_names:
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
				if resname == 'ptdins': resname_actual = ptdins_resname
				else: resname_actual = resname
				anum = analysis_names.index(aname)
				mset = msets[anum]
				#---fit
				msds_mean = mean(msds_master[anum][mset.resnames.index(resname_actual)],axis=0)
				datbuf = [int(len(msds_mean)*.1)+2,int(len(msds_mean)*.9)+1]
				[q,r], covariance = polyfit(time_master[anum][datbuf[0]:datbuf[1]],
					msds_mean[datbuf[0]:datbuf[1]],1,cov=True)
				print 'Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
				#---plot all curves
				if show_means:
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
						label=ion_label+extra_label_list[anum]+'\n'+\
							'D = '+'{0:.2f}'.format(q*10**6/4.)+\
							r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
					elif resname_group == 'protonation':
						if resname_actual != ptdins_resname:
							label = 'with '+ion_label+' and '+\
								extra_label_list[anum]+'\n'+\
									'D = '+'{0:.2f}'.format(q*10**6/4.)+\
									r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
						else:
							label = extra_label_list[anum]+' and '+ion_label+\
								'\n'+'D = '+'{0:.2f}'.format(q*10**6/4.)+\
									r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
					else:
						label=ion_label+', D = '+'{0:.2f}'.format(q*10**6/4.)+\
							r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
					ax.plot(time_master[anum][datbuf[0]:datbuf[1]],msds_mean[datbuf[0]:datbuf[1]],
						color=color,lw=3,
						label=label)
					ax.plot(time_master[anum],msds_mean,'--',color=color,lw=3)
			#---plot settings
			ax.legend(loc='upper left',fontsize=fsaxlegend-2)
			ax.set_xscale('log')
			ax.set_yscale('log')
			if fix_lims != None:
				ax.set_xlim(fix_lims[0])
				ax.set_ylim(fix_lims[1])
			if 0: ax.grid(b=True,which='major',color='k',linestyle='-')
			if 0: ax.grid(b=True,which='minor',color='k',linestyle='--',alpha=0.2)
			ax.grid(True)
			ax.set_xlabel('Time (ps)',fontsize=fsaxlabel)
			plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
			plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper'))
			if pnum < len(valid_resnames)-1:
				ax.set_xticklabels([])
				ax.set_xlabel('')
			if 0: props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
			props = dict(boxstyle='round',facecolor='w',alpha=1.0,lw=1)
			ax.text(0.6,0.8,('PtdIns' if (resname_group == 'protonation' and resname != resname_actual) else 
				proper_residue_names_long[(ptdins_resname if resname != resname_actual else resname)]),
				fontsize=fsaxlabel+2,
				transform=ax.transAxes,horizontalalignment='center',bbox=props)
			if pnum == 0: ax.set_title('Lipid diffusion',fontsize=fsaxtitle)

		#---single y-axis label	
		fig.text(0.035, 0.5,'MSD (nm$^{2}$)',rotation="vertical",va="center",
			fontsize=fsaxtitle,backgroundcolor='w',alpha=1.)
	
		#---save
		plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-lipid-diffusion-SUMMARY-'+\
			'-'.join(analysis_names)+'.png',dpi=300)
		if showplots: plt.show()
		plt.close(fig)
		
if 'plot_lipid_diffusion_ptdins' in routine:

	#---figure
	fig = plt.figure()
	gs = gridspec.GridSpec(len(valid_resnames),1,wspace=0.0,hspace=0.0)
	#---settings
	show_means = True
	fix_lims = [(5*10**2,100000),(10**-2,2*10**0)]
	#---temporary setting
	if resname_group == 'protonation': fix_lims = [(5*10**1,100000),(10**-2,2*10**0)]
	
	valid_resnames = ['ptdins']
	resname = valid_resnames[0]

	pnum = valid_resnames.index(resname)
	ax = fig.add_subplot(gs[pnum])
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		if resname == 'ptdins': resname_actual = ptdins_resname
		else: resname_actual = resname
		anum = analysis_names.index(aname)
		mset = msets[anum]
		#---fit
		msds_mean = mean(msds_master[anum][mset.resnames.index(resname_actual)],axis=0)
		datbuf = [int(len(msds_mean)*.1)+2,int(len(msds_mean)*.9)+1]
		[q,r], covariance = polyfit(time_master[anum][datbuf[0]:datbuf[1]],
			msds_mean[datbuf[0]:datbuf[1]],1,cov=True)
		print 'Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
		#---plot all curves
		if show_means:
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
				label=ion_label+extra_label_list[anum]+'\n'+\
					'D = '+'{0:.2f}'.format(q*10**6/4.)+\
					r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
			elif resname_group == 'protonation':
				if resname_actual != ptdins_resname:
					label = 'with '+ion_label+' and '+\
						extra_label_list[anum]+'\n'+\
							'D = '+'{0:.2f}'.format(q*10**6/4.)+\
							r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
				else:
					label = extra_label_list[anum]+' and '+ion_label+\
						'\n'+'D = '+'{0:.2f}'.format(q*10**6/4.)+\
							r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'		
			else:
				label=ion_label+', D = '+'{0:.2f}'.format(q*10**6/4.)+\
					r'$\:\mathrm{{\mathbf{\mu} m}^2{s}^{-1}}$'
			ax.plot(time_master[anum][datbuf[0]:datbuf[1]],msds_mean[datbuf[0]:datbuf[1]],
				color=color,lw=3,
				label=label)
			ax.plot(time_master[anum],msds_mean,'--',color=color,lw=3)
	#---plot settings
	ax.legend(loc='upper left',fontsize=fsaxlegend-2)
	ax.set_xscale('log')
	ax.set_yscale('log')
	if fix_lims != None:
		ax.set_xlim(fix_lims[0])
		ax.set_ylim(fix_lims[1])
	if 0: ax.grid(b=True,which='major',color='k',linestyle='-')
	if 0: ax.grid(b=True,which='minor',color='k',linestyle='--',alpha=0.2)
	ax.grid(True)
	ax.set_xlabel('Time (ps)',fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper'))
	if pnum < len(valid_resnames)-1:
		ax.set_xticklabels([])
		ax.set_xlabel('')
	if 0: props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
	props = dict(boxstyle='round',facecolor='w',alpha=1.0,lw=1)
	'''
	if resname_group == 'protonation':
		ax.set_title('PtdIns(4,5)P$_2^{-3}$ or  PtdIns(4,5)P$_2^{-4}$ or PtdIns(4,5)P$_2^{-5}$',
			fontsize=fsaxlabel)
	elif resname_group == 'phosphate_position':
		ax.set_title('PtdIns(4,5)P$_2^{-3}$ or  PtdIns(4,5)P$_2^{-4}$ or PtdIns(4,5)P$_2^{-5}$',
			fontsize=fsaxlabel)	
	'''
	if resname_group == 'protonation': ax.set_title('PtdIns Protonation States',fontsize=fsaxlabel)
	elif resname_group == 'phosphate_position': ax.set_title('PtdIns Phosphate Position',fontsize=fsaxlabel)		
	#---single y-axis label	
	ax.set_ylabel('MSD (nm$^{2}$)',fontsize=fsaxlabel)

	#---save
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-lipid-diffusion-single-'+\
		'-'.join(analysis_names)+'.png',dpi=300)
	if showplots: plt.show()
	plt.close(fig)
		

#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

if 'batch_override' not in globals():

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
		][:3]
	
	#---alternate tests
	if compare_phosphate_position:
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]

	#---script functions
	routine = [
		'calc',
		'plot',
		'calc_lipid_ion',
		][-1]
	
	#---specify a single pair for analysis, see batch if you want more
	pairing = [
		['ptdins','ptdins'],
		['POPC','POPC'],
		['ptdins','DOPE'],
		['ptdins','DOPS'],
		['ptdins','DOPC'],
		['ptdins','CHL1'],
		][0]
		
	#---pair specifications for lipid-ion g(r)
	pairings_lipid_ion = [
		['ptdins','ion'],
		][:]
	
	#---set the g(r) bin width in Angstroms
	binsizeabs = 1

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if type(routine) == str: routine = [routine]

if 'calc' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		status('status: pairing = '+str(pairing)+' system = '+aname+'\n')
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			#---note: removed a check to see if the pkl exists from this space
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)

			#---radial distribution calculation
			fast_pbc = True
			pair = pairing
			#---removed loop over pairings and outsourced to a batch script
			result_data = MembraneData('gr2d',label='-'.join(pair))
			allcurves = []
			for mononum in range(2):
				for lnum in range(2):
					if pair[lnum] == 'ptdins': 
						lipid_selector = 'resname '+ptdins_resname+' and '+key_atom_selector[ptdins_resname]
					else:
						lipid_selector = 'resname '+pair[lnum]+' and '+key_atom_selector[pair[lnum]]
					allselect_lipids = mset.universe.selectAtoms(lipid_selector)
					validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
						set([i-1 for i in allselect_lipids.resids()])))
					mset.selections.append(sum([allselect_lipids.residues[
						list(allselect_lipids.resids()).index(i+1)].selectAtoms(lipid_selector) 
						for i in validresids]))
				allcurves_by_monolayer = []
				if whichframes == None: frameslice = range(len(mset.universe.trajectory))
				else: frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
				if mset.selections[0] == 0.0 or mset.selections[1] == 0.0:
					status('status: missing lipids in this monolayer\n')
					allcurves_by_monolayer = [[] for i in range(len(frameslice))]
				else:
					for frameno in frameslice:
						status('status: frame = '+str(frameno+1)+'/'+str(len(frameslice)))
						mset.gotoframe(frameno)
						vecs = mset.vec(frameno)
						pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
						pts2 = array(mset.get_points(frameno,selection_index=1))[:,0:2]
						if fast_pbc:
							pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
							pts2 = array(mset.get_points(frameno,selection_index=1))[:,0:2]
							dmat2a = scipy.spatial.distance.cdist(array([[i] 
								for i in pts1[:,0]]),array([[i] for i in pts2[:,0]]))
							dmat2b = scipy.spatial.distance.cdist(array([[i] 
								for i in pts1[:,1]]),array([[i] for i in pts2[:,1]]))
							dmat2a_pbc = dmat2a-1*(dmat2a>vecs[0]/2.)*vecs[0]
							dmat2b_pbc = dmat2b-1*(dmat2b>vecs[1]/2.)*vecs[1]
							pairdists = sqrt(dmat2a_pbc**2+dmat2b_pbc**2)
							#---note that this fast method has trouble with the square-circle regions
							cutoff = sqrt(sum((vecs[:2]/2.)**2))
							cutoff = mean(vecs[:2])/2.
						else:
							st = time.time()
							status('status: frame = '+str(frameno+1)+'/'+str(len(frameslice)))
							mset.gotoframe(frameno)
							vecs = mset.vec(frameno)
							pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
							pts2 = array(mset.get_points(frameno,selection_index=1))[:,0:2]
							points = pts2
							dims=[0,1]
							ans = []
							for p in points:
								for tr in [[i,j] for i in arange(-2,2+1,1) for j in arange(-2,2+1,1) 
									if not (i == 0 and j == 0)]:
									ans.append([p[i]+tr[i]*vecs[i] 
										if i in dims else p[i] for i in range(2)])
							pts2pbc = concatenate((points,array(ans)))
							pairdists = scipy.spatial.distance.cdist(pts1,pts2pbc)
							pts2 = pts2pbc
							cutoff = 2*(vecs[0] if vecs[0]<vecs[1] else vecs[1])
						tmp = pairdists.min(axis=1)
						sysarea = pi*cutoff**2
						scanrange = arange(0,int(cutoff),binsizeabs)
						histcomb = []
						for r in range(len(pairdists)):
							row = pairdists[r]
							hist,binedge = numpy.histogram(array(row)[array(row)!=0.],
								range=(0.,int(max(scanrange))),bins=len(scanrange)-1)
							mid = (binedge[1:]+binedge[:-1])/2
							histcomb.append(hist)
						grcurve = sum(histcomb,axis=0)/float(len(histcomb))
						allcurves_by_monolayer.append(grcurve)
					cutoff = 2*min([int(i) for i in np.mean(mset.vecs,axis=0)[0:2]])
				allcurves.append(allcurves_by_monolayer)
				mset.selections = []
			if 'pts1' not in globals(): break
			points_counts = [shape(pts1),shape(pts2)]
			pair_selects = pair[0:2]
			pair_name = '-'.join(pair)
			#---load up the membranedata object with calculation details
			savelist = ['cutoff','sysarea','points_counts','binsizeabs','pair_selects','pair_name']
			for item in savelist:
				result_data.addnote([item,globals()[item]])
			result_data.data = [[allcurves[0][i],allcurves[1][i]] for i in range(len(frameslice))]
			result_data.label = pair_name
			mset.store = [result_data]
			pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
			#---modify names if the frame selection doesn't perfectly match the original gmx timeslice
			if '-'.join(aname.split('-')[1:]) != \
				'-'.join(specname_pickle(sysname,trajfile[0]).split('.')[-1:]):
				specname_mod = '.'.join(specname_pickle(sysname,traj).split('.')[:-1])+'.'+\
					'-'.join(aname.split('-')[1:])
			else: specname_mod = specname_pickle(sysname,traj)
			pickledump(mset,'pkl.gr2d.'+'-'.join(pairnames)+'.'+\
				specname_mod+'.pkl',directory=pickles)
			checktime()
			del mset
			del result_data
			
if 'calc_lipid_ion' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		status('status: pairing = '+str(pairing)+' system = '+aname+'\n')
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			#---note: removed a check to see if the pkl exists from this space
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)

			#---load the ions
			mset_ions = MembraneSet()
			grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
				keysysname='ions_sysname',keytrajsel='ions_trajsel')
			mset_ions.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')

			#---radial distribution calculation
			fast_pbc = True
			pair = pairing
			#---removed loop over pairings and outsourced to a batch script
			result_data = MembraneData('gr2d',label='-'.join(pair))
			allcurves = []
			for mononum in range(2):
				for lnum in range(2):
					if pair[lnum] == 'ion':
						allselect_lipids = mset_ions.universe.selectAtoms('name '+str(ion_name))
						mset_ions.selections.append(allselect_lipids)
					else:
						if pair[lnum] == 'ptdins': 
							lipid_selector = 'resname '+ptdins_resname+' and '+key_atom_selector[ptdins_resname]
						elif pair[lnum] == 'ions':
							 lipid_selector = 'name '+str(ion_name)
						else:
							lipid_selector = 'resname '+pair[lnum]+' and '+key_atom_selector[pair[lnum]]
						allselect_lipids = mset.universe.selectAtoms(lipid_selector)
						validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
							set([i-1 for i in allselect_lipids.resids()])))
						mset.selections.append(sum([allselect_lipids.residues[
							list(allselect_lipids.resids()).index(i+1)].selectAtoms(lipid_selector) 
							for i in validresids]))
				allcurves_by_monolayer = []
				if whichframes == None: frameslice = range(len(mset.universe.trajectory))
				else: frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
				if mset.selections[0] == 0.0 or mset_ions.selections[0] == 0.0:
					status('status: missing lipids in this monolayer\n')
					allcurves_by_monolayer = [[] for i in range(len(frameslice))]
				else:
					for frameno in frameslice:
						status('status: frame = '+str(frameno+1)+'/'+str(len(frameslice)))
						mset.gotoframe(frameno)
						mset_ions.gotoframe(frameno)
						vecs = mset.vec(frameno)
						pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
						pts2 = array(mset_ions.get_points(frameno,selection_index=0))[:,0:2]
						if fast_pbc:
							pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
							pts2 = array(mset_ions.get_points(frameno,selection_index=0))[:,0:2]
							dmat2a = scipy.spatial.distance.cdist(array([[i] 
								for i in pts1[:,0]]),array([[i] for i in pts2[:,0]]))
							dmat2b = scipy.spatial.distance.cdist(array([[i] 
								for i in pts1[:,1]]),array([[i] for i in pts2[:,1]]))
							dmat2a_pbc = dmat2a-1*(dmat2a>vecs[0]/2.)*vecs[0]
							dmat2b_pbc = dmat2b-1*(dmat2b>vecs[1]/2.)*vecs[1]
							pairdists = sqrt(dmat2a_pbc**2+dmat2b_pbc**2)
							cutoff = sqrt(sum((vecs[:2]/2.)**2))
							cutoff = mean(vecs[:2])/2.
						tmp = pairdists.min(axis=1)
						sysarea = pi*cutoff**2
						scanrange = arange(0,int(cutoff),binsizeabs)
						histcomb = []
						for r in range(len(pairdists)):
							row = pairdists[r]
							hist,binedge = numpy.histogram(array(row)[array(row)!=0.],
								range=(0.,int(max(scanrange))),bins=len(scanrange)-1)
							mid = (binedge[1:]+binedge[:-1])/2
							histcomb.append(hist)
						grcurve = sum(histcomb,axis=0)/float(len(histcomb))
						allcurves_by_monolayer.append(grcurve)
					cutoff = 2*min([int(i) for i in np.mean(mset.vecs,axis=0)[0:2]])
				allcurves.append(allcurves_by_monolayer)
				mset.selections = []
			if 'pts1' not in globals(): break
			points_counts = [shape(pts1),shape(pts2)]
			pair_selects = pair[0:2]
			pair_name = '-'.join(pair)
			#---load up the membranedata object with calculation details
			savelist = ['cutoff','sysarea','points_counts','binsizeabs','pair_selects','pair_name']
			for item in savelist:
				result_data.addnote([item,globals()[item]])
			result_data.data = [[allcurves[0][i],allcurves[1][i]] for i in range(len(frameslice))]
			result_data.label = pair_name
			mset.store = [result_data]
			pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
			#---modify names if the frame selection doesn't perfectly match the original gmx timeslice
			if '-'.join(aname.split('-')[1:]) != \
				'-'.join(specname_pickle(sysname,trajfile[0]).split('.')[-1:]):
				specname_mod = '.'.join(specname_pickle(sysname,traj).split('.')[:-1])+'.'+\
					'-'.join(aname.split('-')[1:])
			else: specname_mod = specname_pickle(sysname,traj)
			pickledump(mset,'pkl.gr2d.'+'-'.join(pairnames)+'.'+\
				specname_mod+'.pkl',directory=pickles)
			checktime()
			del mset
			del result_data

#---LOADS
#-------------------------------------------------------------------------------------------------------------

if 'calc' not in routine and ('plot' in routine or 'plot_lipid_ion' in routine) and 'msets' not in globals():
	msets = []
	#---only analyze single pairs for now
	pair = pairing
	print 'in calc'
	print pair
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
		#---modify names if the frame selection doesn't perfectly match the original gmx timeslice
		if '-'.join(aname.split('-')[1:]) != \
			'-'.join(specname_pickle(sysname,trajfile[0]).split('.')[-1:]):
			specname_mod = '.'.join(specname_pickle(sysname,trajfile[0]).split('.')[:-1])+'.'+\
				'-'.join(aname.split('-')[1:])
		else: specname_mod = specname_pickle(sysname,trajfile[0])
		mset = unpickle(pickles+'pkl.gr2d.'+'-'.join(pairnames)+'.'+\
			specname_mod+'.pkl')
		msets.append(mset)


#---PLOT
#-------------------------------------------------------------------------------------------------------------

if 'plot' in routine:	
	print 'start '
	print pair
	print pairing
	pair = pairing
	#---check if the comparison and pair are valid (even though empty pkl may exist, this is easier)
	if (not any([mset == None for mset in msets]) and \
		any([all([i in mset.resnames for i in pairnames]) for mset in msets])) or \
		resname_group == 'protonation':

		#---seemingly redundant routine to ensure the right bins
		all_cutoff_bins = []
		for aname in analysis_names:
			print aname
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			mset = msets[analysis_names.index(aname)]
			mdat = mset.store[0]
			pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
			#---handle cases where lipids may only be on one monolayer
			allcurves_both_monolayers = []	
			if not all([i in mset.resnames for i in pairnames]): break
			which_monolayers = [i for i in range(2) if 
				(len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])])>0 and
				len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[1])])>0)]
			for i in which_monolayers:
				allcurves_both_monolayers.append(mdat.get(['monolayer',i]))
			allcurves = concatenate([i for i in allcurves_both_monolayers])
			cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
			all_cutoff_bins.append(cutoff_bins)
		consensus_cutoff_bins = min(all_cutoff_bins)
		maxpeak = 0	
		
		fig = plt.figure(figsize=(8,6))
		gs = mpl.gridspec.GridSpec(1,1)
		ax = fig.add_subplot(gs[0])
		plotted_objects = 0
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
			mset = msets[analysis_names.index(aname)]
			mdat = mset.store[0]
			#---handle cases where lipids may only be on one monolayer
			allcurves_both_monolayers = []	
			if not all([i in mset.resnames for i in pairnames]) and resname_group != 'protonation': break
			which_monolayers = [i for i in range(2) if 
				(len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])])>0 and
				len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[1])])>0)]
			for i in which_monolayers:
				allcurves_both_monolayers.append(mdat.get(['monolayer',i]))
			allcurves = concatenate([i for i in allcurves_both_monolayers])
			binsizeabs = mdat.getnote('binsizeabs')
			cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
			scanrange = arange(0,int(consensus_cutoff_bins))
			nbins = len(scanrange)-1
			avgcurv = np.mean(array([i[0:nbins] for i in allcurves]),axis=0)
			hist,binedge = numpy.histogram([1 for i in range(1000)],
				range=(0,max(scanrange)*binsizeabs),bins=len(scanrange)-1)
			mid = (binedge[1:]+binedge[:-1])/2
			binwidth = binsizeabs
			areas = [pi*binwidth*mid[i]*2 for i in range(len(mid))]
			#---previously used mdat.getnote('points_counts')[0][0]
			nlipids = max([i[0] for i in mdat.getnote('points_counts')])
			vecs = np.mean(mset.vecs,axis=0)
			#---color selection section
			if resname_group == 'phosphate_position':
				color = color_dictionary_aamd(ionname=ion_name,lipid_resname=ptdins_resname,
					comparison='ions_phospate_position')
			elif resname_group == 'protonation':
				color = color_dictionary_aamd(ionname=ion_name,
					lipid_resname=ptdins_resname,comparison='protonation')
			else:
				color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
			if 'color_overrides' in globals(): color = color_overrides[analysis_names.index(aname)]
			if 'marker_overrides' in globals(): plot_marker = marker_overrides[analysis_names.index(aname)]
			else: plot_marker = '-'
			labelname = [proper_residue_names[i] for i in pairnames]
			if pairnames[0] == pairnames[1]:
				label = labelname[0]+' (self), '+ion_label+\
					(early_late[analysis_names.index(aname)] if 'early_late' in globals() else '')+\
					(extra_label_list[analysis_names.index(aname)] 
						if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
			else:
				label = labelname[0]+'-'+labelname[1]+', '+ion_label+\
					(early_late[analysis_names.index(aname)] if 'early_late' in globals() else '')+\
					(extra_label_list[analysis_names.index(aname)] 
						if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
			#---plot
			if len(avgcurv) > 0: plotted_objects += 1
			gr_curve = avgcurv/areas/(nlipids/(vecs[0]*vecs[1]))
			maxpeak = max(gr_curve) if maxpeak < max(gr_curve) else maxpeak
			ax.plot(mid/10.,gr_curve,plot_marker,lw=2.5,color=color,
				label=label)
		ax.set_ylim((0,(2.0 if maxpeak > 1.5 else 1.5)))
		ax.legend()
		ax.grid(True)
		ax.set_xlabel(r'$\mathbf{r}\:\mathrm{(nm)}$', fontsize=fsaxlabel)
		ax.set_ylabel(r'$\mathrm{g(\mathbf{r})}$', fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxticks)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxticks)
		ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
		plt.legend(loc='lower right')
		print plotted_objects
		print pairnames
		#---save
		if plotted_objects > 0 :
			plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-gr2d-'+\
				'-'.join(pairnames)+'.'+\
				'-'.join(analysis_names)+'.png',\
				dpi=300,bbox_inches='tight')
		if showplots: plt.show()
		plt.close(fig)

if 'plot_lipid_ion' in routine:	

	#---seemingly redundant routine to ensure the right bins
	all_cutoff_bins = []
	for aname in analysis_names:
		print aname
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		mset = msets[analysis_names.index(aname)]
		mdat = mset.store[0]
		print pairnames
		pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
		#---handle cases where lipids may only be on one monolayer
		allcurves_both_monolayers = []	
		which_monolayers = [i for i in range(2) if 
			(len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])])>0)]
		for i in which_monolayers:
			allcurves_both_monolayers.append(mdat.get(['monolayer',i]))
		allcurves = concatenate([i for i in allcurves_both_monolayers])
		cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
		all_cutoff_bins.append(cutoff_bins)
	consensus_cutoff_bins = min(all_cutoff_bins)
	maxpeak = 0	

	fig = plt.figure(figsize=(8,6))
	gs = mpl.gridspec.GridSpec(1,1)
	ax = fig.add_subplot(gs[0])
	plotted_objects = 0
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
		mset = msets[analysis_names.index(aname)]
		mdat = mset.store[0]
		#---handle cases where lipids may only be on one monolayer
		allcurves_both_monolayers = []	
		which_monolayers = [i for i in range(2) if 
			(len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])])>0)]
		for i in which_monolayers:
			allcurves_both_monolayers.append(mdat.get(['monolayer',i]))
		allcurves = concatenate([i for i in allcurves_both_monolayers])
		binsizeabs = mdat.getnote('binsizeabs')
		cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
		scanrange = arange(0,int(consensus_cutoff_bins))
		nbins = len(scanrange)-1
		avgcurv = np.mean(array([i[0:nbins] for i in allcurves]),axis=0)
		hist,binedge = numpy.histogram([1 for i in range(1000)],
			range=(0,max(scanrange)*binsizeabs),bins=len(scanrange)-1)
		mid = (binedge[1:]+binedge[:-1])/2
		binwidth = binsizeabs
		areas = [pi*binwidth*mid[i]*2 for i in range(len(mid))]
		#---previously used mdat.getnote('points_counts')[0][0]
		nlipids = max([i[0] for i in mdat.getnote('points_counts')])
		vecs = np.mean(mset.vecs,axis=0)
		#---color selection section
		if resname_group == 'phosphate_position':
			color = color_dictionary_aamd(ionname=ion_name,lipid_resname=ptdins_resname,
				comparison='ions_phospate_position')
		elif resname_group == 'protonation':
			color = color_dictionary_aamd(ionname=ion_name,
				lipid_resname=ptdins_resname,comparison='protonation')
		else:
			color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
		if 'color_overrides' in globals(): color = color_overrides[analysis_names.index(aname)]
		if 'marker_overrides' in globals(): plot_marker = marker_overrides[analysis_names.index(aname)]
		else: plot_marker = '-'
		labelname = [(ion_label if i == 'ion' else proper_residue_names[i]) for i in pairnames]
		if pairnames[0] == pairnames[1]:
			label = labelname[0]+' (self)'+\
				(early_late[analysis_names.index(aname)] if 'early_late' in globals() else '')+\
				(extra_label_list[analysis_names.index(aname)] 
					if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
		else:
			label = labelname[0]+', '+labelname[1]+\
				(early_late[analysis_names.index(aname)] if 'early_late' in globals() else '')+\
				(extra_label_list[analysis_names.index(aname)] 
					if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
		#---plot
		if len(avgcurv) > 0: plotted_objects += 1
		gr_curve = avgcurv/areas/(nlipids/(vecs[0]*vecs[1]))
		maxpeak = max(gr_curve) if maxpeak < max(gr_curve) else maxpeak
		ax.plot(mid/10.,gr_curve,plot_marker,lw=2.5,color=color,
			label=label)
	ax.set_ylim((0,(3.0 if maxpeak > 1.5 else 1.5)))
	ax.set_xlim((0,4.))
	ax.legend()
	ax.grid(True)
	ax.set_xlabel(r'$\mathbf{r}\:\mathrm{(nm)}$', fontsize=fsaxlabel)
	ax.set_ylabel(r'$\mathrm{g(\mathbf{r})}$', fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxticks)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxticks)
	ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
	plt.legend(loc='upper right')
	#---save
	if plotted_objects > 0 :
		plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-gr2d-'+\
			'-'.join(pairnames)+'.'+\
			'-'.join(analysis_names)+'.png',\
			dpi=300,bbox_inches='tight')
	if showplots: plt.show()
	plt.close(fig)


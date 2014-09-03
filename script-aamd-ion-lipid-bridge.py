#!/usr/bin/python -i

from membrainrunner import *

execfile('locations.py')
execfile('header-aamd.py')

# TODO:
# - Save to pickle
# - Keep track of which oxygens the ion is bridging
# - Add angle dependency
# - Add this to a batch run script

if 'batch_override' not in globals():
	# Testing non-batch version on v532, PtdIns-Cal
	analysis_names = [
		                 'v530-40000-90000-50',
		                 'v531-40000-90000-50-bridge',
		                 'v532-40000-90000-50-bridge',
		                 'v509-40000-90000-50',
		                 'v510-40000-90000-50',
		                 'v511-40000-90000-50',
		                 'v533-40000-90000-50',
		                 'v534-40000-90000-50',
	                 ][1:3]

	routine = [
		          'calculate',
		          'plot'
	          ][0:1]

	# This is not going to change; this is all script-aamd-ion-lipid-bridge can do.
	pairings_lipid_ion = [
		                     ['ptdins', 'ion'],
	                     ][:]

	# Binding distance in Angstrom, angle dependence not coded yet.
	binding_cutoff = 4.0

	fast_pbc = True  # I'm not even sure what fast_pbc = False does

if type(routine) == str: routine = [routine]

if 'calculate' in routine:

	for aname in analysis_names:
		status('status: Pairing = ' + str(pairings_lipid_ion[0]) + ' system = ' + aname + '\n')

		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile, trajfile = trajectory_lookup(analysis_descriptors, aname, globals())

		for traj in trajfile:
			# Load lipids and split into monolayers
			mset = MembraneSet()
			mset.load_trajectory((basedir + '/' + grofile, basedir + '/' + traj), resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			# Load ions
			mset_ions = MembraneSet()
			# grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			# keysysname='ions_sysname',keytrajsel='ions_trajsel')
			# Definition of trajectory_lookup in membrainrunner has the order of the arguments reversed
			grofile, trajfile = trajectory_lookup(analysis_descriptors, aname, globals(),
			                                      keytrajsel='ions_trajsel', keysysname='ions_sysname')
			mset_ions.load_trajectory((basedir + '/' + grofile, basedir + '/' + trajfile[0]), resolution='aamd')
			# Unwrap ions
			fast_pbc = True
			# This is probably going to be overidden in a batch script, so hack it for now:
			pair = pairings_lipid_ion[0]
			result_data = MembraneData('bridge', label='-'.join(pair))
			num_ions_binding = []  # Binding to at least one lipid
			num_ions_bridging = []  # Binding to at least two lipids
			num_lipids_ion_binding = []  # This array holds a list of the number of lipids all ions are binding
			# by frame. The range is set by variable k.
			which_oxygens = []
			which_o_pairs = []
			total_oxygens = 0

			# Loop over the lipid and the ion selections
			for lnum in range(2):
				if pair[lnum] == 'ion':
					# Choosing only cation for analysis, set by ion_name in header-aamd.py dictionary
					# Heads up: this nomenclature is confusing because allselect_lipids is used for Cal.
					allselect_lipids = mset_ions.universe.selectAtoms('name ' + str(ion_name))
					# Check v532, Cal: 308 atoms, CL: 396 atoms
					mset_ions.selections.append(allselect_lipids)
				else:
					# For this analysis:
					# a) restrict to PtdIns
					# b) use phosphate oxygen atoms, not phosphate phosphorus
					# Make a new entry in header-aamd.py dictionary for this
					# keylipidatoms becomes binding_atoms in the dictionary
					# To use non-PtdIns, need to make new time-slices.
					if pair[lnum] == 'ptdins':
						lipid_selector = 'resname ' + ptdins_resname + ' and ' + binding_atoms[ptdins_resname]
					elif pair[lnum] == 'ions':
						lipid_selector = 'name ' + str(ion_name)
					else:
						lipid_selector = 'resname ' + pair[lnum] + ' and ' + binding_atoms[pair[lnum]]
					phosphodiester_check = list(mset.universe.selectAtoms('resname ' + ptdins_resname + \
					                                                      ' and (name O14 or name O13)'))
					if len(phosphodiester_check) == 0:
						status('status: Atoms O13 and O14 are not in the structure file. You are probably " \
						"using a structure file from /keylipidatoms/ which does not have them included. You can " \
						"only check for binding to inositol phosphate groups and not the phosphodiester.')

					allselect_lipids = mset.universe.selectAtoms(lipid_selector)
					# Instead of doing this loop residue-by-residue, this will tell us how to group the array pts1
					num_lipids = mset.universe.selectAtoms(lipid_selector).numberOfResidues()
					#validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
					#	set([i-1 for i in allselect_lipids.resids()])))
					#mset.selections.append(sum([allselect_lipids.residues[
					#	list(allselect_lipids.resids()).index(i+1)].selectAtoms(lipid_selector)
					#	for i in validresids]))
					mset.selections.append(allselect_lipids)
			# Select the frames for analysis (default = all)
			if whichframes == None:
				frameslice = range(len(mset.universe.trajectory))
			else:
				frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
			for frameno in frameslice:
				status('status: frame = ' + str(frameno + 1) + '/' + str(len(frameslice)))
				mset.gotoframe(frameno)
				mset_ions.gotoframe(frameno)
				vecs = mset.vec(frameno)
				# Not sure what selection_index is doing in the get_points() function
				# mset.selections has two elements, not sure what second entry (0.0) is:
				# In [19]: mset.selections
				# Out[19]: [<AtomGroup with 40 atoms>, 0.0]
				pts1 = array(mset.get_points(frameno, selection_index=0))
				pts2 = array(mset_ions.get_points(frameno, selection_index=0))
				# All the binding oxygen coordinates in current frame
				# There are len(pts1)/num_lipids oxygens per lipid
				# So when checking for ion distances to two different lipids, make sure they are
				# more than len(pts1)/num_lipids apart.
				pts1 = array(mset.get_points(frameno, selection_index=0))
				# All the ion coordinates in current frame
				pts2 = array(mset_ions.get_points(frameno, selection_index=0))
				# Could do pre-filtering if I weren't lazy.
				distance_x = scipy.spatial.distance.cdist(array([[i]
				                                                 for i in pts1[:, 0]]),
				                                          array([[i] for i in pts2[:, 0]]))
				distance_y = scipy.spatial.distance.cdist(array([[i]
				                                                 for i in pts1[:, 1]]),
				                                          array([[i] for i in pts2[:, 1]]))
				distance_z = scipy.spatial.distance.cdist(array([[i]
				                                                 for i in pts1[:, 2]]),
				                                          array([[i] for i in pts2[:, 2]]))
				unwrapped_x = distance_x - 1 * (distance_x > vecs[0] / 2.) * vecs[0]
				unwrapped_y = distance_y - 1 * (distance_y > vecs[1] / 2.) * vecs[1]
				unwrapped_z = distance_z - 1 * (distance_z > vecs[2] / 2.) * vecs[2]
				distance_matrix = sqrt(unwrapped_x ** 2 + unwrapped_y ** 2 + unwrapped_z ** 2)
				# Since we know the number of oxygens per lipid for this analysis, we can now
				# *quickly* group the distances by lipid instead of writing a wrapper loop to do
				# by residue analysis
				points_per_lipid = len(pts1) / num_lipids
				num_binding_sites = len(pts1)
				num_ions = len(pts2)
				# i is ion number (e.g 308), j is atom number for all atoms in all lipids (e.g. 240),
				ion_to_lipid = [[distance_matrix[:, i][j:j + points_per_lipid] \
				                 for j in arange(0, num_binding_sites, points_per_lipid)] \
				                for i in range(num_ions)]

				# min_ion_to_lipid[ion][:] shows the minimum distance between a given ion and any binding
				# oxygen on a given lipid
				min_ion_to_lipid = [[ion_to_lipid[i][j][:].min() \
				                     for j in range(num_lipids)] \
				                    for i in range(num_ions)]

				# If we keep track of *where* the minimum is, then we can also keep track of which oxygen
				# is doing the binding. And if we keep track of both minima, then we know which two oxygens
				# are being bridged by the ion. This assumes two factors
				# a) That SelectAtoms() outputs coordinates in the same order as inputs
				# b) That cdist stores the distance between input i and input j in the ij'th entry in output matrix

				min_index_and_value = [[min(enumerate(ion_to_lipid[i][j]), key=operator.itemgetter(1)) \
				                        for j in range(num_lipids)] \
				                       for i in range(num_ions)]


				# Check if a given ion is binding to a single lipid:
				binding_truth_table = [any([i < binding_cutoff for i in min_ion_to_lipid[j]]) \
				                       for j in range(num_ions)]
				fraction_binding = float(sum(binding_truth_table)) / float(len(binding_truth_table))
				bridging_truth_table = [len(where([i < binding_cutoff for i in min_ion_to_lipid[j]])[0]) >= 2 \
				                        for j in range(num_ions)]
				fraction_bridging = float(sum(bridging_truth_table)) / float(len(bridging_truth_table))

				# Need to normalize this by the number of PI2P that /could/ be bridged, by finding
				# the closest distance between oxygens of two adjacent residues and comparing with binding_cuytoff.

				num_ions_binding.append(sum(binding_truth_table))
				num_ions_bridging.append(sum(bridging_truth_table))
				this_frame = []  # This array holds a list of the number of lipids an ion is binding /per frame/
				which_oxygens_frame = []  # This array holds a list of which oxygens an ion is binding /per frame/
				which_o_pairs_frame = []  # This array holds a list of oxygen pairs an ion is binding /per frame/

				# Assumption: no ion is coordinating more than 4 lipids
				for k in range(4):
					# The first method works just fine, but does not use the matrix that keeps track of the oxygens
					# num_lipids_binding = [len(where([i < binding_cutoff for i in min_ion_to_lipid[j]])[0]) == k \
					#                      for j in range(num_ions)]
					ion_to_lipid_min_dists = [[row[1] for row in min_index_and_value[j]] for j in range(num_ions)]
					num_lipids_binding = [len(where([i < binding_cutoff for i in ion_to_lipid_min_dists[j]])[0])
					                      == k for j in range(num_ions)]
					# This array holds the index of the oxygen on a specific lipid (0 to 6 or 8) that is closest to a
					# given ion
					ion_to_lipid_index_of_min = [[row[0] for row in min_index_and_value[0]] for j in range(num_ions)]
					if k == 1:
						# print 'There should be '+str(len(where(num_lipids_binding)[0]))+' single oxygen-ion ' \
						#	   'interactions in this frame.'
						total_oxygens += len(where(num_lipids_binding)[0])
						which_ions_binding_single_lipid = where(num_lipids_binding)  # Where it is true
						for i in which_ions_binding_single_lipid[0]:
							this_ion = ion_to_lipid_min_dists[i]
							# Find location of minimum
							min_index = this_ion.index(min(this_ion))
							this_oxygen = ion_to_lipid_index_of_min[i][min_index]
							which_oxygens_frame.append(this_oxygen)
					elif k == 2:
						# Here we want to keep track of *pairs* not just single atoms.
						total_oxygens += len(where(num_lipids_binding)[0])
						which_ions_binding_two_lipids = where(num_lipids_binding)
						for i in which_ions_binding_two_lipids[0]:
							this_ion = ion_to_lipid_min_dists[i]
							min_indices = argsort(ion_to_lipid_min_dists[i])[0:2]
							these_oxygens = ion_to_lipid_index_of_min[i][min_indices[0]], ion_to_lipid_index_of_min[i][
								min_indices[1]]
							which_o_pairs_frame.append(these_oxygens)
					this_frame.append(sum(num_lipids_binding))
				num_lipids_ion_binding.append(this_frame)
				which_oxygens.append(which_oxygens_frame)
				which_o_pairs.append(which_o_pairs_frame)

			oxygen_pairs = [item for sublist in which_o_pairs for item in sublist]
			from collections import Counter
			pair_sorted = []
			for pair in oxygen_pairs:
				# Sorting the entries should get rid of the double (x,y) and (y,x) counts
				pair_sorted.append(tuple(sort(pair)))
			count = Counter()
			for pair in pair_sorted:
				count[pair] += 1


			# Using a dictionary would be nice, but D.keys() returns items in a random unsorted way.
			# D = {'OP52':0, 'OP53':1, 'OP54':2, 'OP42':3,'OP43':4, 'OP44':5, 'O13':6, 'O14':7}
			D = ['OP52', 'OP53', 'OP54', 'OP42', 'OP43', 'OP44', 'O13', 'O14']
			Descriptors = ['P5', 'P5', 'P5', 'P4', 'P4', 'P4', 'P1', 'P1']
			interactive = 1
		if interactive == 1:
			checktime()
			print 'Counting oxygens'
			big_list = [item for sublist in which_oxygens for item in sublist]
			# This really ought to be replaced because it's so slow...
			oxygen_count = dict((i, big_list.count(i)) for i in big_list)
			plt.bar([i for i in range(len(oxygen_count))],
			        [float(oxygen_count[i]) / len(big_list) for i in range(len(oxygen_count))], align='center')
			plt.xticks(range(len(D)), [D[i] for i in range(len(D))])
			plt.title('Distribution of oxygens binding to '+str(ion_name)+' ions which only bind one oxygen')
			plt.show()

			# Can do some more data de-duplication by collecting all P5/P4/P1 pairs together.

			most_common = count.most_common()[0:10]
			plt.bar([x for x in range(len(most_common))], [float(most_common[i][1])/sum(count.values()) for i in range(len(most_common))])
			plt.xticks([x for x in range(len(most_common))], \
			           [D[most_common[i][0][0]]+str('-')+D[most_common[i][0][1]] for i in range(len(most_common))] )
			#plt.xticks([x for x in range(len(most_common))], \
			#           [Descriptors[most_common[i][0][0]]+str('-')+Descriptors[most_common[i][0][1]] for i in range(len(most_common))] )

			plt.xticks(rotation=45)
			plt.ylabel('Fraction of all '+str(ion_name)+' bridges')
			plt.title('Distribution of oxygen pairs coordinated by '+str(ion_name)+ ' ions binding two lipids')
			plt.show()

		# The following bits help figure out what's going wrong when the results are donkey bonkers...
		debug = 0
		if debug:
			import itertools as it

			d = [[max([y - x for x, y in it.combinations(ion_to_lipid[i][j][:], 2)]) \
			      for j in range(num_lipids)] \
			     for i in range(num_ions)]
			print "The maximum difference between the distance of all oxygen atoms of a specific lipid " \
			      "to given ion can be checked by accessing d[ion][lipid] and this difference should " \
			      "be relatively small (probably less than 10 A)."
			# For each ion, how far away is nearest lipid binding oxygen?
			# distance_matrix.min(axis=0)
			# For each ion, how far away are the /two/ nearest lipid binding oxygens?
			# distance_matrix[distance_matrix[:,i].argsort()[0:2],i]
			# For each lipid binding oxygen, how far away is nearest ion?
			# distance_matrix.min(axis=1)
			# Minimum distance from ion to any lipid binding oxygen:
			distance_to_lipid = [distance_matrix[:, i].min() for i in range(shape(distance_matrix)[1])]
			# Check if a given ion is binding to /any/ oxygen atom:
			[distance_to_lipid[i] < binding_cutoff for i in range(len(distance_to_lipid))]

	# Store the parameters and the results in membraindata and a pickle
	points_counts = [shape(pts1), shape(pts2)]
	pair_selects = pair[0:2]
	pair_name = '-'.join(pair)
	savelist = ['binding_cutoff', 'points_counts', 'pair_selects', 'pair_name']
	for item in savelist:
		result_data.addnote([item, globals()[item]])
	result_data.data = [num_lipids_ion_binding]
	result_data.label = pair_name
	mset.store = [result_data]
	pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
	# ---modify names if the frame selection doesn't perfectly match the original gmx timeslice
	if '-'.join(aname.split('-')[1:]) != '-'.join(specname_pickle(sysname, trajfile[0]).split('.')[-1:]):
		specname_mod = '.'.join(specname_pickle(sysname, traj).split('.')[:-1]) + '.' + '-'.join(aname.split('-')[1:])
	else:
		specname_mod = specname_pickle(sysname, traj)
	pickledump(mset.store[0], 'pkl.bridge.' + '-'.join(pairnames) + '.' + specname_mod + '.pkl', directory=pickles)
	checktime()
	del mset
	del mset_ions
	del result_data

if 'plot' in routine:
	pair = pairings_lipid_ion[0]
	pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in pair]
	# ---check if the comparison and pair are valid (even though empty pkl may exist, this is easier)
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
			                    (len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])]) > 0 and
			                     len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[1])]) > 0)]
			for i in which_monolayers:
				allcurves_both_monolayers.append(mdat.get(['monolayer', i]))
			allcurves = concatenate([i for i in allcurves_both_monolayers])
			cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
			all_cutoff_bins.append(cutoff_bins)
		consensus_cutoff_bins = min(all_cutoff_bins)
		maxpeak = 0
		fig = plt.figure(figsize=(8, 6))
		gs = mpl.gridspec.GridSpec(1, 1)
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
			                    (len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[0])]) > 0 and
			                     len(mset.monolayer_by_resid[i][mset.resnames.index(pairnames[1])]) > 0)]
			for i in which_monolayers:
				allcurves_both_monolayers.append(mdat.get(['monolayer', i]))
			allcurves = concatenate([i for i in allcurves_both_monolayers])
			binsizeabs = mdat.getnote('binsizeabs')
			cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
			scanrange = arange(0, int(consensus_cutoff_bins))
			nbins = len(scanrange) - 1
			avgcurv = np.mean(array([i[0:nbins] for i in allcurves]), axis=0)
			hist, binedge = numpy.histogram([1 for i in range(1000)],
			                                range=(0, max(scanrange) * binsizeabs), bins=len(scanrange) - 1)
			mid = (binedge[1:] + binedge[:-1]) / 2
			binwidth = binsizeabs
			areas = [pi * binwidth * mid[i] * 2 for i in range(len(mid))]
			#---previously used mdat.getnote('points_counts')[0][0]
			nlipids = max([i[0] for i in mdat.getnote('points_counts')])
			vecs = np.mean(mset.vecs, axis=0)
			#---color selection section
			if resname_group == 'phosphate_position':
				color = color_dictionary_aamd(ionname=ion_name, lipid_resname=ptdins_resname,
				                              comparison='ions_phospate_position')
			elif resname_group == 'protonation':
				color = color_dictionary_aamd(ionname=ion_name,
				                              lipid_resname=ptdins_resname, comparison='protonation')
			else:
				color = color_dictionary_aamd(ionname=ion_name, comparison='ions')
			if 'color_overrides' in globals(): color = color_overrides[analysis_names.index(aname)]
			if 'marker_overrides' in globals():
				plot_marker = marker_overrides[analysis_names.index(aname)]
			else:
				plot_marker = '-'
			labelname = [proper_residue_names[i] for i in pairnames]
			if pairnames[0] == pairnames[1]:
				label = labelname[0] + ' (self), ' + ion_label + \
				        (early_late[analysis_names.index(aname)] if 'early_late' in globals() else '') + \
				        (extra_label_list[analysis_names.index(aname)]
				         if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
			else:
				label = labelname[0] + '-' + labelname[1] + ', ' + ion_label + \
				        (early_late[analysis_names.index(aname)] if 'early_late' in globals() else '') + \
				        (extra_label_list[analysis_names.index(aname)]
				         if ('ptdins' not in pair and 'extra_label_list' in globals()) else '')
			#---plot
			if len(avgcurv) > 0: plotted_objects += 1
			gr_curve = avgcurv / areas / (nlipids / (vecs[0] * vecs[1]))
			maxpeak = max(gr_curve) if maxpeak < max(gr_curve) else maxpeak
			ax.plot(mid / 10., gr_curve, plot_marker, lw=2.5, color=color,
			        label=label)
		ax.set_ylim((0, (2.0 if maxpeak > 1.5 else 1.5)))
		ax.legend()
		ax.grid(True)
		ax.set_xlabel(r'$\mathbf{r}\:\mathrm{(nm)}$', fontsize=fsaxlabel)
		ax.set_ylabel(r'$\mathrm{g(\mathbf{r})}$', fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(), fontsize=fsaxticks)
		plt.setp(ax.get_yticklabels(), fontsize=fsaxticks)
		ax.set_title(composition_name + ' bilayer', fontsize=fsaxtitle)
		plt.legend(loc='lower right')
		print plotted_objects
		print pairnames
		#---save
		if plotted_objects > 0:
			plt.savefig(pickles + '/PREPARE-FIGURES/' + 'fig-gr2d-' + \
			            '-'.join(pairnames) + '.' + \
			            '-'.join(analysis_names) + '.png', \
			            dpi=300, bbox_inches='tight')
		if showplots: plt.show()
		plt.close(fig)

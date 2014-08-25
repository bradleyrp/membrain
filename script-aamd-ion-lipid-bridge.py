#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

if 'batch_override' not in globals():

	compare_phosphate_position = False
	if compare_phosphate_position:
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	# Testing non-batch version on v532, PtdIns-Cal
	analysis_names = [
		'v530-40000-90000-50',
		'v531-40000-90000-50',
		'v532-40000-90000-50-bridge',
		'v509-40000-90000-50',
		'v510-40000-90000-50',
		'v511-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		][2:3]

	routine = [
		'calc_lipid_ion',
		][0:1]

	# This is not going to change
	pairings_lipid_ion = [
		['ptdins','ion'],
		][:]

	# Size in Angstroms for binning position data
	binsizeabs = 1

	# Binding distance in Angstrom, angle dependence not coded yet.
	binding_cutoff = 3.0


if type(routine) == str: routine = [routine]

if 'calc_lipid_ion' in routine:

	for aname in analysis_names:
		status('status: Pairing = '+str(pairings_lipid_ion[0])+' system = '+aname+'\n')

		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())

		for traj in trajfile:
			# Load lipids and split into monolayers
			mset = MembraneSet()
			mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
			mset.identify_monolayers(director)
			if residues == 'infer':
				residues = list(set(mset.universe.residues.resnames()))
			mset.identify_residues(residues)
			# Load ions
			mset_ions = MembraneSet()
			# grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			#	keysysname='ions_sysname',keytrajsel='ions_trajsel')
			# Definition of trajectory_lookup in membrainrunner has the order of the arguments reversed
			grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
				keytrajsel='ions_trajsel',keysysname='ions_sysname')
			mset_ions.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
			# Unwrap ions
			fast_pbc = True
			# This is probably going to be overidden in a batch script, so hack it for now:
			pair = pairings_lipid_ion[0]
			# I'll have to create a new data type and figure out whether to put it in membrainrunner or membraindata...
			# result_data = MembraneData('bridge',label='-'.join(pair))
			num_ions_binding = []  # Binding to at least one lipid
			num_ions_bridging = [] # Binding to at least two lipids
			num_lipids_ion_binding = [] # This array simply holds a list of the number of lipids a given ion is binding
			# Loop over the lipid and the ion selections
			for lnum in range(2):
				if pair[lnum] == 'ion':
					# Choosing only cation for analysis, set by ion_name in header-aamd.py dictionary
					# Heads up: this nomenclature is confusing because allselect_lipids is used for Cal.
					allselect_lipids = mset_ions.universe.selectAtoms('name '+str(ion_name))
					# Check v532, Cal: 308 atoms, CL: 396 atoms
					mset_ions.selections.append(allselect_lipids)
				else:
					# For this analysis:
					# a) restrict to PtdIns
					# b) use phosphate oxygen atoms, not phosphate phosphorus
					# Make a new entry in header-aamd.py dictionary for this
					# keylipidatoms becomes binding_atoms in the dictionary
					if pair[lnum] == 'ptdins':
						lipid_selector = 'resname '+ptdins_resname+' and '+binding_atoms[ptdins_resname]
					elif pair[lnum] == 'ions':
						lipid_selector = 'name '+str(ion_name)
					else:
						lipid_selector = 'resname '+pair[lnum]+' and '+binding_atoms[pair[lnum]]
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
			# allcurves_by_monolayer = []
			# Select the frames for analysis (default = all)
			if whichframes == None: frameslice = range(len(mset.universe.trajectory))
			else: frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
			for frameno in frameslice:
				status('status: frame = '+str(frameno+1)+'/'+str(len(frameslice)))
				mset.gotoframe(frameno)
				mset_ions.gotoframe(frameno)
				vecs = mset.vec(frameno)
				# Not sure what selection_index is doing in the get_points() function
				# mset.selections has two elements, not sure what second entry (0.0) is:
				# In [19]: mset.selections
				# Out[19]: [<AtomGroup with 40 atoms>, 0.0]
				pts1 = array(mset.get_points(frameno,selection_index=0))
				pts2 = array(mset_ions.get_points(frameno,selection_index=0))
				if fast_pbc: # Not sure what fast_pbc=False does...
					# All the binding oxygen coordinates in current frame
					# There are len(pts1)/num_lipids oxygens per lipid
					# So when checking for ion distances to two different lipids, make sure they are
					# more than len(pts1)/num_lipids apart.
					pts1 = array(mset.get_points(frameno,selection_index=0))
					# All the ion coordinates in current frame
					pts2 = array(mset_ions.get_points(frameno,selection_index=0))
					# Could do pre-filtering if I weren't lazy.
					distance_x = scipy.spatial.distance.cdist(array([[i]
						for i in pts1[:,0]]),array([[i] for i in pts2[:,0]]))
					distance_y = scipy.spatial.distance.cdist(array([[i]
						for i in pts1[:,1]]),array([[i] for i in pts2[:,1]]))
					distance_z = scipy.spatial.distance.cdist(array([[i]
						for i in pts1[:,2]]),array([[i] for i in pts2[:,2]]))
					unwrapped_x = distance_x-1*(distance_x>vecs[0]/2.)*vecs[0]
					unwrapped_y = distance_y-1*(distance_y>vecs[1]/2.)*vecs[1]
					unwrapped_z = distance_z-1*(distance_z>vecs[2]/2.)*vecs[2]
					distance_matrix = sqrt(unwrapped_x**2+unwrapped_y**2+unwrapped_z**2)
					# Since we know the number of oxygens per lipid for this analysis, we can now
					# *quickly* group the distances by lipid instead of writing a wrapper loop to do
					# by residue analysis
					points_per_lipid = len(pts1)/num_lipids
					num_binding_sites = len(pts1)
					num_ions = len(pts2)
					# i is ion number (e.g 308), j is atom number for all atoms in all lipids (e.g. 240),
					ion_to_lipid = [[distance_matrix[:,i][j:j+points_per_lipid] \
							for j in arange(0,num_binding_sites,points_per_lipid)] \
					        for i in range(num_ions)]
					# min_ion_to_lipid[ion][:] shows the minimum distance between a given ion and any binding
					# oxygen on a given lipid
					min_ion_to_lipid = [[ion_to_lipid[i][j][:].min() \
		                    for j in range(num_lipids)] \
                            for i in range(num_ions)]
					# Check if a given ion is binding to a single lipid:
					binding_truth_table = [any([i<binding_cutoff for i in min_ion_to_lipid[j]]) \
							for j in range(num_ions)]
					fraction_binding = float(sum(binding_truth_table))/float(len(binding_truth_table))

					bridging_truth_table = [len(where([i<binding_cutoff for i in min_ion_to_lipid[j]])[0]) >= 2 \
					                  for j in range(num_ions)]
					fraction_bridging = float(sum(bridging_truth_table))/float(len(bridging_truth_table))

					num_ions_binding.append(sum(binding_truth_table))
					num_ions_bridging.append(sum(bridging_truth_table))
					this_frame = [] # This array holds a list of the number of lipids an ion is binding /per frame/
					for k in range(4):
						num_lipids_binding = [len(where([i < binding_cutoff for i in min_ion_to_lipid[j]])[0]) == k \
						                      for j in range(num_ions)]
						# num_lipids_ion_binding.append([k,sum(num_lipids_binding)])
						binding = [k, sum(num_lipids_binding)]
						this_frame.append(binding)
				num_lipids_ion_binding.append(this_frame)
			# Sanity check:
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
				distance_to_lipid = [distance_matrix[:,i].min() for i in range(shape(distance_matrix)[1])]
				# Check if a given ion is binding to /any/ oxygen atom:
				[distance_to_lipid[i]<binding_cutoff for i in range(len(distance_to_lipid))]

			# Keep this in mind for later:
			cutoff = sqrt(sum((vecs[:]/2.)**2))
			cutoff = mean(vecs[:])/2.
			system_area = pi*cutoff**2
			scan_range = arange(0,int(cutoff),binsizeabs)

if 'plot' in routine:
	for aname in analysis_names:
		status('status: Plotting = '+str(pairings_lipid_ion[0])+' system = '+aname+'\n')


'''
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
'''

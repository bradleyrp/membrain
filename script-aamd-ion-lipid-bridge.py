#!/usr/bin/python -i

from membrainrunner import *
from collections import Counter


execfile('locations.py')
execfile('header-aamd.py')

# TODO:
# - Add angle dependency
# - Split up plot/compute, save pickle, work on plot code and batch

if 'batch_override' not in globals():
	# Testing non-batch version on v532, PtdIns-Cal
	analysis_names = [
		                 'v530-40000-90000-50',
		                 'v531-40000-90000-50-bridge',
		                 'v532-40000-90000-50-bridge',
		                 'v509-40000-90000-50-bridge',
		                 'v510-40000-90000-50-bridge',
		                 'v511-40000-90000-50-bridge',
		                 'v533-40000-90000-50',
		                 'v534-40000-90000-50',
	                 ][1:6]

	routine = [
		          'calculate',
		          'plot'
	          ][1:2]

	# This is not going to change; this is all script-aamd-ion-lipid-bridge can do.
	pairings_lipid_ion = [
		                     ['ptdins', 'ion'],
	                     ][:]

	# Binding distance in Angstrom, angle dependence not coded yet.
	# We have to find a good compromise between Cal and Mag
	binding_cutoff = 4.0

	fast_pbc = True  # I'm not even sure what fast_pbc = False does

	# Using a dictionary would be nice, but D.keys() returns items in a random unsorted way.
	# D = {'OP52':0, 'OP53':1, 'OP54':2, 'OP42':3,'OP43':4, 'OP44':5, 'O13':6, 'O14':7}
	D = ['OP52', 'OP53', 'OP54', 'OP42', 'OP43', 'OP44', 'O13', 'O14']
	Descriptors = ['P5', 'P5', 'P5', 'P4', 'P4', 'P4', 'P1', 'P1']


def pbc_unwrap(pts1, pts2):
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
	
	return distance_matrix

def get_lipid_ion_coords (frameno):
	pts1 = array(mset.get_points(frameno, selection_index=0))
	pts2 = array(mset_ions.get_points(frameno, selection_index=0))
	distance_matrix = pbc_unwrap(pts1,pts2)
	points_per_lipid = len(pts1) / num_lipids
	num_binding_sites = len(pts1)

	return (pts1, pts2, points_per_lipid, distance_matrix, num_binding_sites)

def distance_to_lipids (distance_matrix, points_per_lipid, num_binding_sites):
	# i is ion number (e.g 308), j is atom number for all binding sites in all lipids (e.g. 240),
	ion_to_lipid = [[distance_matrix[:, i][j:j + points_per_lipid] \
	                 for j in arange(0, num_binding_sites, points_per_lipid)] \
	                for i in range(num_ions)]

	# min_ion_to_lipid[ion][:] shows the minimum distance between a given ion and any binding site
	# on a given lipid
	# In [9]: shape(min_ion_to_lipid)
	# Out[9]: (308, 40) = (num_ions, num_lipids)

	min_ion_to_lipid = [[ion_to_lipid[i][j][:].min() \
	                     for j in range(num_lipids)] \
	                    for i in range(num_ions)]

	# If we keep track of *where* the minimum is, then we can also keep track of which oxygen
	# is doing the binding. And if we keep track of both minima, then we know which two oxygens
	# are being bridged by the ion. This assumes two factors:
	# a) That SelectAtoms() outputs coordinates in the same order as inputs
	# b) That cdist stores the distance between input i and input j in the ij'th entry in output matrix

	min_index_and_value = [[min(enumerate(ion_to_lipid[i][j]), key=operator.itemgetter(1)) \
	                        for j in range(num_lipids)] \
	                       for i in range(num_ions)]

	return (ion_to_lipid, min_ion_to_lipid, min_index_and_value)


def load(structure, gro, traj, lipids):
	structure.load_trajectory((basedir + '/' + gro, basedir + '/' + traj), resolution='aamd')
	if lipids:
		structure.identify_monolayers(director)
		residues = list(set(structure.universe.residues.resnames()))
		structure.identify_residues(residues)
		lipid_selector = 'resname ' + ptdins_resname + ' and ' + binding_atoms[ptdins_resname]
		phosphodiester_check = list(structure.universe.selectAtoms('resname ' + ptdins_resname + \
	                                                      ' and (name O14 or name O13)'))
		if len(phosphodiester_check) == 0:
			status('status: Atoms O13 and O14 are not in the structure file. You are probably " \
			"using a structure file from /keylipidatoms/ which does not have them included. You can " \
			"only check for binding to inositol phosphate groups and not the phosphodiester.')
		selection = structure.universe.selectAtoms(lipid_selector)
		num = structure.universe.selectAtoms(lipid_selector).numberOfResidues()
	else:
		ion_selector = 'name ' + str(ion_name)
		selection = structure.universe.selectAtoms(ion_selector)
		num = structure.universe.selectAtoms(ion_selector).numberOfResidues()
	structure.selections.append(selection)
	return (structure, num)

def define_binding ():
	if str(ion_name) == 'MG':
		binding_cutoff = 4.0
	elif str(ion_name) == 'Cal':
		binding_cutoff = 3.0
	else:
		binding_cutoff = 4.0
		print 'Check sodium binding range.'
	return(binding_cutoff)

def binding_number ():
	# Stuff
	print 'Not implemented'

def save_data (result_data, oxygens, oxygen_count, oxygen_pairs, pair_sorted, count, pairs_within_cutoff):
	pair_name = ptdins_resname +'-' + ion_name
	savelist = ['pair_name', 'binding_cutoff', 'D'] # These are metadata.
	for item in savelist:
		result_data.addnote([item, globals()[item]])
	result_data.data = [oxygens, oxygen_count, oxygen_pairs, pair_sorted, count, pairs_within_cutoff]
	result_data.label = ['Index of oxygens binding to ion per frame, see notes for corresponding Dictionary', \
	                     'Counter() object tracking index of oxygens binding to ion', \
	                     'Index of oxygen pairs binding to ion per frame', \
	                     'Sorted index of oxygen pairs binding to ion', \
                         'Counter() object tracking index of oxygen pairs binding to ion', \
	                     'Array containing the number of pairs within binding_cutoff for each frame'
	]
	mset.store = [result_data]
	pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in analysis_pair]
	# ---modify names if the frame selection doesn't perfectly match the original gmx timeslice
	if '-'.join(aname.split('-')[1:]) != '-'.join(specname_pickle(sysname, lipid_trajfile[0]).split('.')[-1:]):
		specname_mod = '.'.join(specname_pickle(sysname, lipid_trajfile[0]).split('.')[:-1]) + '.' + '-'.join(aname.split('-')[1:])
	else:
		specname_mod = specname_pickle(sysname, lipid_trajfile[0])
	pickledump(mset.store[0], 'pkl.bridge.' + '-'.join(pairnames) + '.' + specname_mod + '.pkl', directory=pickles)

if type(routine) == str: routine = [routine]

if 'calculate' in routine:

	for aname in analysis_names:
		status('status: Pairing = ' + str(pairings_lipid_ion[0]) + ' system = ' + aname + '\n')

		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		lipid_grofile, lipid_trajfile = trajectory_lookup(analysis_descriptors, aname, globals())
		ion_grofile, ion_trajfile =  trajectory_lookup(analysis_descriptors, aname, globals(),
			                                      keytrajsel='ions_trajsel', keysysname='ions_sysname')
		for traj in lipid_trajfile:
			mset, num_lipids = load(MembraneSet(), lipid_grofile, lipid_trajfile[0], lipids=True)
			mset_ions, num_ions = load(MembraneSet(), ion_grofile, ion_trajfile[0], lipids=False)
			binding_cutoff = define_binding()

			# This is probably going to be overidden in a batch script, so hack it for now:
			analysis_pair = pairings_lipid_ion[0]

			# Containers
			###########################################################
			result_data = MembraneData('bridge', label='-'.join(analysis_pair))
			num_ions_binding = []  # Binding to at least one lipid
			num_ions_bridging = []  # Binding to at least two lipids
			num_lipids_ion_binding = []  # This array holds a list of the number of lipids all ions are binding
			# by frame. The range is set by variable k.
			which_oxygens = [] # Array containing which oxygens bind to given ion
			which_o_pairs = [] # Array containing which two oxygens are bridged by ion
			pairs_within_cutoff = [] # Array containing the number of residues within binding_cutoff
			total_oxygens = 0 # Error checking code to make sure we catch all interacting oxygens
			############################################################

			# Select the frames for analysis (default = all)
			if whichframes == None:
				frameslice = range(len(mset.universe.trajectory))
			else:
				frameslice = [i for i in range(len(mset.universe.trajectory)) if i in whichframes]
			for frameno in frameslice:
			#for frameno in range(999,1000):
				status('status: frame = ' + str(frameno + 1) + '/' + str(len(frameslice)))
				mset.gotoframe(frameno)
				mset_ions.gotoframe(frameno)
				vecs = mset.vec(frameno)

				pts1, pts2, points_per_lipid, distance_matrix, num_binding_sites = get_lipid_ion_coords(frameno)
				ion_to_lipid, min_ion_to_lipid, min_index_and_value = distance_to_lipids(distance_matrix, points_per_lipid, num_binding_sites)

				# Check if a given ion is binding to a single lipid:
				binding_truth_table = [any([i < binding_cutoff for i in min_ion_to_lipid[j]]) \
				                       for j in range(num_ions)]
				fraction_binding = float(sum(binding_truth_table)) / float(len(binding_truth_table))
				bridging_truth_table = [len(where([i < binding_cutoff for i in min_ion_to_lipid[j]])[0]) >= 2 \
				                        for j in range(num_ions)]
				fraction_bridging = float(sum(bridging_truth_table)) / float(len(bridging_truth_table))

				# Need to normalize this by the number of PI2P that /could/ be bridged, by finding
				# the closest distance between oxygens of two adjacent residues and comparing with binding_cutoff.

				num_ions_binding.append(sum(binding_truth_table))
				num_ions_bridging.append(sum(bridging_truth_table))
				this_frame = []  # This array holds a list of the number of lipids an ion is binding /per frame/
				which_oxygens_frame = []  # This array holds a list of which oxygens an ion is binding /per frame/
				which_o_pairs_frame = []  # This array holds a list of oxygen pairs an ion is binding /per frame/
				pairs_within_cutoff_frame = set() # This set holds how many pairs (40*39/2) are within binding cutoff /per frame/

				# I think this is already calculated:
				# ion_to_lipid_min_dists = [[row[1] for row in min_index_and_value[j]] for j in range(num_ions)]
				ion_to_lipid_min_dists = min_ion_to_lipid
				# Assumption: no ion is coordinating more than 3 lipids
				for k in range(4):
					num_lipids_binding = [len(where([i < binding_cutoff for i in min_ion_to_lipid[j]])[0])
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
							#print 'Ion '+str(i)+' is bridging residues '+str(min_indices)
							these_oxygens = ion_to_lipid_index_of_min[i][min_indices[0]], \
							                ion_to_lipid_index_of_min[i][min_indices[1]]
							which_o_pairs_frame.append(these_oxygens)
					this_frame.append(sum(num_lipids_binding))

				#################################################################
				# Calculate minimum lipid-lipid distances for normalization
				# There *must* be a much faster way to do this, but I can't figure it out right now.
				# One speedup would be not wasting memory on squareform, but then it's more challenging
				# to check if two elements come from the same residue.
				# We really only need the upper triangular part, sigh.
				# lipid_coords = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(pts1))
				lipid_coords = pbc_unwrap(pts1, pts1)
				for i in range(len(lipid_coords)):
					for j in range(len(lipid_coords)):
						if lipid_coords[i,j] < 2*binding_cutoff:
							if i / points_per_lipid != j / points_per_lipid: # This integer division should check if
								# they are coming from the same residue.
									if tuple(sort([i/points_per_lipid,j/points_per_lipid])) not in pairs_within_cutoff_frame:
										# I am so wary this is working...
										pairs_within_cutoff_frame.add(tuple(sort([i/points_per_lipid,j/points_per_lipid])))
										# from mayavi import mlab
										# mlab.points3d(pts1[:,0],pts1[:,1],pts1[:,2],scale_factor=0.7)
										#print 'Residue '+str(i/points_per_lipid)+' is close to residue '+str(j/points_per_lipid)
				#################################################################

				pairs_within_cutoff.append(len(pairs_within_cutoff_frame))
				num_lipids_ion_binding.append(this_frame)
				which_oxygens.append(which_oxygens_frame)
				which_o_pairs.append(which_o_pairs_frame)

			#################################################################
			# Post-processing each trajectory
			# Single lipid "binding"
			oxygens = [item for sublist in which_oxygens for item in sublist]
			oxygen_count = Counter()
			for oxygen in oxygens:
				oxygen_count[oxygen] += 1
			# Two lipid "bridging"
			oxygen_pairs = [item for sublist in which_o_pairs for item in sublist]
			pair_sorted = []
			for pair in oxygen_pairs:
				# Sorting the entries should get rid of separately counting (x,y) and (y,x)
				pair_sorted.append(tuple(sort(pair)))
			count = Counter()
			for pair in pair_sorted:
				count[pair] += 1
			most_common = count.most_common()[0:10]
			#################################################################


		checktime()
		save_data (result_data, oxygens, oxygen_count, oxygen_pairs, pair_sorted, count, pairs_within_cutoff)
		checktime()

		del mset
		del mset_ions
		del result_data

if 'plot' in routine:
	if 'compute' not in routine:
			analysis_pair = pairings_lipid_ion[0]
			for aname in analysis_names:
				status('status: Plotting= ' + str(pairings_lipid_ion[0]) + ' system = ' + aname + '\n')
				# This pairname piece of junk is buggy as fsck.
				# The following block is not strictly necessary, but is not buggy.
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
				lipid_grofile, lipid_trajfile = trajectory_lookup(analysis_descriptors, aname, globals())
				ion_grofile, ion_trajfile =  trajectory_lookup(analysis_descriptors, aname, globals(),
					                                      keytrajsel='ions_trajsel', keysysname='ions_sysname')
				for traj in lipid_trajfile:
					mset, num_lipids = load(MembraneSet(), lipid_grofile, lipid_trajfile[0], lipids=True)
					mset_ions, num_ions = load(MembraneSet(), ion_grofile, ion_trajfile[0], lipids=False)
					binding_cutoff = define_binding()
					analysis_pair = pairings_lipid_ion[0]
				pairnames = [(i if i != 'ptdins' else ptdins_resname) for i in analysis_pair]
				# ---modify names if the frame selection doesn't perfectly match the original gmx timeslice
				if '-'.join(aname.split('-')[1:]) != '-'.join(specname_pickle(sysname, lipid_trajfile[0]).split('.')[-1:]):
					specname_mod = '.'.join(specname_pickle(sysname, lipid_trajfile[0]).split('.')[:-1]) + '.' + '-'.join(aname.split('-')[1:])
				else:
					specname_mod = specname_pickle(sysname, lipid_trajfile[0])
				results = unpickle(pickles+'pkl.bridge.' + '-'.join(pairnames) + '.' + specname_mod + '.pkl')

				(oxygens, oxygen_count, oxygen_pairs, pair_sorted, count, pairs_within_cutoff) = \
					(results.data[0],results.data[1],results.data[2],results.data[3],results.data[4],results.data[5])

				fig = plt.figure(figsize=(8, 6))
				gs = mpl.gridspec.GridSpec(1, 1)
				ax = fig.add_subplot(gs[0])

				plt.bar([i for i in range(len(oxygen_count))],
					        [float(oxygen_count[i]) / len(oxygen_count) for i in range(len(oxygen_count))], align='center')
				plt.xticks(range(len(D)), [D[i] for i in range(len(D))])
				plt.title('Oxygens binding to '+str(ion_name)+' ions binding exactly 1 residue (cutoff = '+str(binding_cutoff)+')' )
				if show: plt.show()

				plt.bar([x for x in range(len(most_common))], [float(most_common[i][1])/sum(count.values()) for i in range(len(most_common))])
				plt.xticks([x for x in range(len(most_common))], \
					           [D[most_common[i][0][0]]+str('-')+D[most_common[i][0][1]] for i in range(len(most_common))] )
				plt.xticks(rotation=45)
				plt.ylabel('Fraction of all '+str(ion_name)+' bridges')
				plt.title('Oxygen pairs coordinated by '+str(ion_name)+ ' ions binding exactly 2 residues (cutoff ='+str(binding_cutoff)+')')


				plt.savefig(pickles + 'fig-bridge.' + '-'.join(pairnames) + '.' + specname_mod + '.png', \
	            dpi=300, bbox_inches='tight')
				if show: plt.show()
				plt.close(fig)
				
	elif 'compute' in routine:
		checktime()
		status('Plotting most recent calculation.')
		show = True
		fig = plt.figure(figsize=(8, 6))
		gs = mpl.gridspec.GridSpec(1, 1)
		ax = fig.add_subplot(gs[0])

		plt.bar([i for i in range(len(oxygen_count))],
			        [float(oxygen_count[i]) / len(oxygen_count) for i in range(len(oxygen_count))], align='center')
		plt.xticks(range(len(D)), [D[i] for i in range(len(D))])
		plt.title('Oxygens binding to '+str(ion_name)+' ions binding exactly 1 residue (cutoff = '+str(binding_cutoff)+')' )
		if show: plt.show()

		plt.bar([x for x in range(len(most_common))], [float(most_common[i][1])/sum(count.values()) for i in range(len(most_common))])
		plt.xticks([x for x in range(len(most_common))], \
			           [D[most_common[i][0][0]]+str('-')+D[most_common[i][0][1]] for i in range(len(most_common))] )
		plt.xticks(rotation=45)
		plt.ylabel('Fraction of all '+str(ion_name)+' bridges')
		plt.title('Oxygen pairs coordinated by '+str(ion_name)+ ' ions binding exactly 2 residues (cutoff ='+str(binding_cutoff)+')')
		if show: plt.show()

		# Plot fraction of all potential bridges occupied over time, histogrammed...

	else:
		print 'No data to plot.'


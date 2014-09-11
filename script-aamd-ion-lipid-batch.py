#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

batch_override = True

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

#---standard selection
analysis_names = [
	'v531-40000-90000-50',
	'v532-40000-90000-50',
	][:]
	
# This is the *only* pairing ion-lipid-bridge can handle now.
pairings_lipid_ion = [
	['ptdins','ion'],
	][:]

all_batches = [
	'calculate',
	'plot',
	][:]

# Binding distance in Angstrom, angle dependence not coded yet.
binding_cutoff = 3.0


#---BATCH
#-------------------------------------------------------------------------------------------------------------

for batch in all_batches:
	if batch == 'calculate':
		routine = 'calculate'
	elif batch == 'calc_lipid_ion':
		routine = 'calc_lipid_ion'		
	elif batch == 'phosphate_position':
		resname_group = 'phosphate_position'
		routine = 'plot' if not lipid_ion else 'plot_lipid_ion'
		compare_phosphate_position = True
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
		extra_label_list = [
			' (with 4,5)',
			' (with 4,5)',
			' (with 3,5)',
			' (with 3,5)',
			]
	elif batch == 'protonation':
		routine = 'plot' if not lipid_ion else 'plot_lipid_ion'
		resname_group = 'protonation'
		compare_phosphate_position = False
		analysis_names = [
			'v509-40000-90000-50',
			'v514-22000-32000-10',
			'v515-20000-30000-10',
			]
		extra_label_list = [
			r' $\mathrm{PIP_2^{-4}}$',
			r' $\mathrm{PIP_2^{-5}}$',
			r' $\mathrm{PIP_2^{-3}}$',
			]
	elif batch == 'phosphate_position_35_only':
		routine = 'plot'
		compare_phosphate_position = True
		analysis_names = [
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
		extra_label_list = [
			' (with 3,5)',
			' (with 3,5)',
			]
	elif batch == 'symmetric':
		routine = 'plot' if not lipid_ion else 'plot_lipid_ion'
		compare_phosphate_position = False
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			]
	elif batch == 'asymmetric':
		routine = 'plot' if not lipid_ion else 'plot_lipid_ion'
		compare_phosphate_position = False
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'early_late':
		routine = 'plot'
		compare_phosphate_position = False
		analysis_names = [
			'v530-40000-65000-50',
			'v530-65000-90000-50',
			]
		early_late = [
			' (early)',
			' (late)',
			]
		color_overrides = ['#708090','k']
		marker_overrides = ['-','-']
	elif batch == 'early_late_v532':
		routine = 'plot'
		compare_phosphate_position = False
		analysis_names = [
			'v532-40000-65000-50',
			'v532-65000-90000-50',
			]
		early_late = [
			' (early)',
			' (late)',
			]
		color_overrides = ['#708090','k']
		marker_overrides = ['-','-']
	else: raise Exception('except: not in batchlist')
	if 'compare_phosphate_position' not in globals(): compare_phosphate_position = False
	if 'resname_group' not in globals(): resname_group = None
	status('status: routine = '+routine+'\n')
	status('status: analysis_names = '+str(analysis_names)+'\n')

	#---loop over requested analyses
	if lipid_ion:
		for pairing in pairings_lipid_ion:
			status('status: pairing = '+str(pairing)+'\n')
			execfile('script-aamd-gr2d.py')
			if 'msets' in globals(): del msets
	else:
		for pairing in pairings:
			status('status: pairing = '+str(pairing)+'\n')
			execfile('script-aamd-gr2d.py')
			if 'msets' in globals(): del msets
	checktime()
	if 'msets' in globals(): del msets
	if 'mset' in globals(): del mset
	if 'extra_label_list' in globals(): del extra_label_list
	if 'early_late' in globals(): del early_late
	if 'marker_overrides' in globals(): del marker_overrides
	if 'color_overrides' in globals(): del color_overrides
	if 'resname_group' in globals(): del resname_group
	

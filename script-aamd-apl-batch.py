#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

batch_override = True

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

#---standard selection
analysis_names = [
	'v530-40000-90000-50',
	'v531-40000-90000-50',
	'v532-40000-90000-50',
	'v509-40000-90000-50',
	'v510-40000-90000-50',
	'v511-40000-90000-50',
	'v514-22000-32000-10',
	'v515-20000-30000-10',
	'v533-40000-90000-50',
	'v534-40000-90000-50',
	'v530-40000-65000-50',
	'v530-65000-90000-50',
	'v534-20000-25000-10',
	][-1:]
	
#---plot settings
nbins = 40
showplots = False
hist_area_max = 140

all_batches = [
	'calc_apl',
	'calc_span',
	'snapshot',
	'phosphate_position_apl',
	'phosphate_position_apl_nochl',
	'asymmetric_apl',
	'asymmetric_apl_nochl',
	'symmetric_apl',
	'phosphate_position_span',
	'asymmetric_span',
	'symmetric_span',
	'protonation_span',
	'summary',
	'protonation',
	'plot_span_angle',
	'plot_span_angle_summary',
	][3:4]

#---BATCH
#-------------------------------------------------------------------------------------------------------------

for batch in all_batches:
	if batch ==	'calc_apl':
		routine = 'calc_apl'
	elif batch ==	'calc_span':
		routine = 'calc_span'
	elif batch == 'phosphate_position_apl':
		routine = 'plot_apl'
		resname_group = 'phosphate_position'
		selector_override = None
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-20000-25000-10',
			'v534-40000-90000-50',
			][:-1]
	elif batch == 'phosphate_position_apl_nochl':
		selector_override = 'no_chl'
		routine = 'plot_apl'
		resname_group = 'phosphate_position'
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	elif batch == 'asymmetric_apl':
		routine = 'plot_apl'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'asymmetric_apl_nochl':
		routine = 'plot_apl'
		selector_override = 'no_chl'
		resname_group = 'all'
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'symmetric_apl':
		routine = 'plot_apl'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			]
	elif batch == 'phosphate_position_span':
		selector_override = None
		routine = 'plot_span'
		resname_group = 'phosphate_position'
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	elif batch == 'asymmetric_span':
		routine = 'plot_span'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'symmetric_span':
		routine = 'plot_span'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			]
	elif batch == 'snapshot':
		routine = 'plot_snapshot'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v530-pbcmol-40000-90000-50',
			]
	elif batch == 'summary':
		routine = 'plot_apl_all_lipids'
		resname_group = 'all'
		selector_override = None
		analysis_names = [
			'v530-40000-90000-50',
			]
	elif batch == 'protonation':
		routine = 'plot_apl'
		resname_group = 'protonation'
		selector_override = None
		analysis_names = [
			'v509-40000-90000-50',
			'v514-22000-32000-10',
			'v515-20000-30000-10',
			]
	elif batch == 'protonation_span':
		routine = 'plot_span'
		resname_group = 'protonation'
		selector_override = None
		analysis_names = [
			'v509-40000-90000-50',
			'v514-22000-32000-10',
			'v515-20000-30000-10',
			]
	elif batch == 'plot_span_angle':
		routine = 'plot_span_angle'
		resname_group = 'all'
		selector_override = None
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
			'v530-40000-65000-50',
			'v530-65000-90000-50',
			]
	elif batch == 'plot_span_angle_summary':
		routine = 'plot_span_angle_summary'
		resname_group = 'all'
		selector_override = None
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
			]
		span_angle_summary_layout = [
			['v530-40000-90000-50','v531-40000-90000-50','v532-40000-90000-50'],
			['v509-40000-90000-50','v510-40000-90000-50','v511-40000-90000-50'],
			['v514-22000-32000-10','v515-20000-30000-10'],
			['v533-40000-90000-50','v534-40000-90000-50'],
			]
		rightlabels = [(analysis_descriptors[aname])['composition_name']+' bilayer' 
			for aname in [
			'v530-40000-90000-50',
			'v509-40000-90000-50',
			'v514-22000-32000-10',
			'v533-40000-90000-50',
			]]
	else: raise Exception('except: not in batchlist')

	flaglist = []
	if 'compare_phosphate_position' not in globals(): compare_phosphate_position = False
	if 'selector_override' not in globals(): selector_override = None
	if selector_override == 'no_chl':
		for aname in analysis_names:
			(analysis_descriptors[aname])['residues'] = residues_aamd_asymmetric_no_chl	
			(analysis_descriptors[aname])['director'] = director_aamd_asymmetric_no_chl
			(analysis_descriptors[aname])['selector'] = selector_aamd_asymmetric_no_chl
		flaglist.append('no_chl')

	status('status: routine = '+routine+'\n')
	status('status: analysis_names = '+str(analysis_names)+'\n')

	#---loop over requested analyses
	execfile('script-aamd-apl.py')
	if 'msets' in globals(): del msets
	checktime()
	

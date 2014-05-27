#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

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
	'v533-40000-90000-50',
	'v534-40000-90000-50',
	'v514-22000-32000-10',
	'v515-20000-30000-10',
	'v530-40000-65000-50',
	'v530-65000-90000-50',
	'v532-40000-65000-50',
	'v532-65000-90000-50',
	][-4:]
	
all_batches = [
	'asymmetric',
	'symmetric',
	'compare_phosphate_position',
	'compare_protonation',
	][-2:-1]

#---BATCH
#-------------------------------------------------------------------------------------------------------------

for batch in all_batches:
	if batch == 'asymmetric':
		routine = [
			'load',
			'compute',
			'plot',
			'plot_lipid_diffusion_summary',
			'plot_lipid_diffusion_detail',
			]
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'symmetric':
		routine = [
			'load',
			'compute',
			'plot',
			'plot_lipid_diffusion_summary',
			'plot_lipid_diffusion_detail',
			]
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			]
	elif batch == 'compare_phosphate_position':
		resname_group = 'phosphate_position'
		routine = [
			'load',
			'compute',
			'plot',
			'plot_lipid_diffusion_summary',
			'plot_lipid_diffusion_detail',
			]
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
	elif batch == 'compare_protonation':
		resname_group = 'protonation'
		routine = [
			'load',
			'compute',
			'plot',
			'plot_lipid_diffusion_summary',
			'plot_lipid_diffusion_detail',
			]
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
	else: raise Exception('except: not in batchlist')
	status('status: analysis_names = '+str(analysis_names)+'\n')
	if 'resname_group' not in globals(): resname_group = None
	execfile('script-aamd-diffusion-lipid.py')
	#---caution: if you comment this line for development, you must uncomment it before running many batches
	#if 'msets' in globals(): del msets
	#if 'resname_group' in globals(): del resname_group
	#if 'extra_label_list' in globals(): del extra_label_list
	checktime()


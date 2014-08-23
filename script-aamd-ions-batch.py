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
	][:1]
	
#---interpolation size 
rounder = 4

#---normalize to bulk concentratio
norm_z_concentration = False

#---develop mode
devmode = True

#---binwidth for ion distributions
binwidth = 1

#---BATCH
#-------------------------------------------------------------------------------------------------------------

for batch in all_batches:
	if batch == 'asymmetric':
		get_ion_alt = True
		routine = [
			'load',
			'compute_z',
			]
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			]
	elif batch == 'symmetric':
		get_ion_alt = True
		routine = [
			'load',
			'compute_z',
			]
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			]
	elif batch == 'compare_phosphate_position':
		resname_group = 'phosphate_position'
		get_ion_alt = True
		routine = [
			'load',
			'compute_z',
			]
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
		extra_label_list = [
			'4,5',
			'4,5',
			'3,5',
			'3,5',
			]
	elif batch == 'compare_protonation':
		resname_group = 'protonation'
		get_ion_alt = True
		routine = [
			'load',
			'compute_z',
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
	status('status: analysis_names = '+str(analysis_names))
	if 'resname_group' not in globals(): resname_group = None
	execfile('script-aamd-ions.py')
	#---caution: if you comment this line for development, you must uncomment it before running many batches
	if not devmode and 'msets' in globals(): del msets
	if not devmode and 'resname_group' in globals(): del resname_group
	if not devmode and 'extra_label_list' in globals(): del extra_label_list
	checktime()


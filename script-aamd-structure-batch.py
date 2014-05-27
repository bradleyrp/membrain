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
	][:]
	
all_batches = [
	'symmetric',
	'asymmetric',
	'phosphate_position',
	'protonation',
	][:]
	
#---settings
showplots = False

#---parameters
rounder = 4

#---peristalsis
peristalsis = False

#---BATCH
#-------------------------------------------------------------------------------------------------------------

for batch in all_batches:
	if batch ==	'asymmetric':
		routine = 'plot'
		compare_phosphate_position = False
		analysis_names = [
			'v530-40000-90000-50',
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			][:]
	elif batch == 'symmetric':
		routine = 'plot'
		compare_phosphate_position = False
		analysis_names = [
			'v509-40000-90000-50',
			'v510-40000-90000-50',
			'v511-40000-90000-50',
			][:]
	elif batch == 'protonation':
		compare_phosphate_position = False
		routine = 'plot'
		analysis_names = [
			'v509-40000-90000-50',
			'v514-22000-32000-10',
			'v515-20000-30000-10',
			][:]
	elif batch == 'phosphate_position':
		compare_phosphate_position = True
		routine = 'plot'
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	else: raise Exception('except: not in batchlist')
	if 'compare_phosphate_position' not in globals(): compare_phosphate_position = False
	status('status: routine = '+routine+'\n')
	status('status: analysis_names = '+str(analysis_names)+'\n')

	#---loop over requested analyses
	execfile('script-aamd-structure.py')
	if 'msets' in globals(): del msets
	checktime()
	

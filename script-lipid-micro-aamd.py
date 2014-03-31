#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

#---plan
analysis_plan = slice(None,None)
analysis_descriptors = {
	'v530-30000-100000-100':
		{'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-director',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.director.xtc',
		'whichframes':None}}
analysis_names = [
	'v530-30000-100000-100'][:]
routine = ['voronoi','tilts','surfnorms'][-1:]


#!/usr/bin/python

'''
deprecated. formerly used in script-lipid-micro.py
'''

from membrainrunner import *
execfile('locations.py')

#---plan
analysis_plan = slice(None,None)
analysis_descriptors = {
	'v612-75000-175000-200':
		{'sysname':'membrane-v612',
		'sysname_lookup':None,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1\,(v2)}$',
		'nprots':1,
		'protein_pkl':None,
		'whichframes':slice(None,None),
		'director':director_cgmd,
		'selector':selector_cgmd,
		'residues':'infer',
		'protein_select':cgmd_protein,
		'testlist':testlist_mono}}
analysis_names = [
	'v612-75000-175000-200'][:]
routine = ['voronoi','tilts','surfnorms'][:1]

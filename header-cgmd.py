#!/usr/bin/python

#---header file containing dictionaries for most CGMD analysis
#---RPB 2014.04.15

if 'analysis_descriptors_extra' not in globals(): analysis_descriptors_extra = {}

#---selections
sel_cgmd_surfacer = ['name PO4 or name POG','name C2A']
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {
	'v614-120000-220000-200': {
		'sysname':'membrane-v614','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'label_text':'4xENTH',
		'nprots':4,
		},
	'v614-40000-140000-200': {
		'sysname':'membrane-v614','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s6-sim-lonestar/md.part0002.40000-140000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'label_text':'4xENTH',
		'nprots':4,
		},
	'v700-500000-600000-200': {
		'sysname':'membrane-v700','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u1-lonestar-longrun/md.part0009.500000-700000-200.xtc',
		'label':r'$\mathrm{{Exo70}\ensuremath{\times}2{\small(para)}}$',
		'timeslice':[500000,600000,200],
		'nprots':2,
		'label_text':'2xExo70(parallel)',
		},
	'v701-60000-160000-200': {
		'sysname':'membrane-v701','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc',
		'label':r'$\mathrm{{Exo70}\ensuremath{\times}2{\small(anti)}}$',
		'nprots':2,
		'label_text':'2xExo70(antiparallel)',
		},
	'v612-75000-175000-200': {
		'sysname':'membrane-v612','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$',
		'label_text':'1xENTH',
		'nprots':1,
		},
	'v612-10000-80000-200': {
		'sysname':'membrane-v612','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-trestles/md.part0003.10000-80000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$',
		'label_text':'1xENTH',
		'nprots':1,
		},
	'v550-400000-500000-160': {
		'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc',
		'label':r'$\mathrm{control}$',
		'label_text':'control',
		'nprots':0,
		},
	'v550-300000-400000-200': {
		'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'label':r'$\mathrm{control}$',
		'label_text':'control',
		'nprots':0,
		},
	'v616-210000-310000-200': {
		'sysname':'membrane-v616','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u6-trestles/md.part0005.210000-310000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'label_text':'8xENTH(close)',
		'nprots':8,
		},
	'v616-110000-209900-200': {
		'sysname':'membrane-v616','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u6-trestles/md.part0004.110000-209900-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'label_text':'8xENTH(close)',
		'nprots':8,
		},
	}
	
#---choose default midplane resolution
gridspacing = 1.0
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'
	
#---coallate two dictionaries
master_dict = dict()
for key in analysis_descriptors.keys():
	if key in analysis_descriptors_extra.keys():
		master_dict.update({key:dict(analysis_descriptors[key],
			**analysis_descriptors_extra[key])})
	else: master_dict.update({key:analysis_descriptors[key]})
analysis_descriptors = master_dict

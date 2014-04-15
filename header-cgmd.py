
#---header file containing dictionaries for most CGMD analysis
#---RPB 2014.04.15

#---selections
sel_cgmd_surfacer = ['name PO4 or name POG','name C2A']
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc'},
	'v614-40000-140000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s6-sim-lonestar/md.part0002.40000-140000-200.xtc'},
	'v700-500000-600000-200':
		{'sysname':'membrane-v700','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u1-lonestar-longrun/md.part0009.500000-700000-200.xtc',
		'timeslice':[500000,600000,200]},
	'v701-60000-160000-200':
		{'sysname':'membrane-v701','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc'},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc'},
	'v612-10000-80000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'s9-trestles/md.part0003.10000-80000-200.xtc'},
	'v550-4000000-500000-160':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc'},
	'v550-300000-400000-200':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc'},
	'v616-61100-117100-1000':
		{'sysname':'membrane-v616','sysname_lookup':None,
		'director':director_cgmd,'selector':selector_cgmd,'protein_select':cgmd_protein,
		'trajsel':'u6-trestles/md.part0004.61100-117100-1000.xtc'},
		}	

#!/usr/bin/python

#---header file containing dictionaries for the DS thesis-round of AAMD analysis
#---RPB 2014.05.20

if 'analysis_descriptors_extra' not in globals(): analysis_descriptors_extra = {}

#---selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
director_aamd_asymmetric_no_chl = ['(name P and not resname CHL1)','(name C218 and not resname CHL1)']
selector_aamd_symmetric = 'name P'
selector_aamd_asymmetric = '(name P and not resname CHL1)'
selector_aamd_asymmetric = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
selector_aamd_asymmetric_no_chl = '(name P and not resname CHL1)'
residues_aamd_symmetric = ['DOPC','DOPS','PI2P']
residues_aamd_asymmetric = ['DOPC','DOPS','DOPE','POPC','PI2P', 'CHL1', 'P35P']
residues_aamd_asymmetric_no_chl = ['DOPC','DOPS','DOPE','POPC','PI2P', 'P35P']
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {

	'v533-40000-54000-100':{
		'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-kraken-md.part0015.40000-54000-100.atomP.xtc',
		'whichframes':None,
		},
	'v534-40000-60000-100':{
		'sysname':'membrane-v534',
		'sysname_lookup':'membrane-v534-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-kraken-md.part0013.40000-60000-100.atomP.xtc',
		'whichframes':None,
		},
	'v530-30000-100000-100':{
		'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.atomP.xtc',
		'whichframes':None,
		},
	'v531-20000-62000-100':{
		'sysname':'membrane-v531',
		'sysname_lookup':'membrane-v531-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-trestles-md.part0007.20000-62000-100.atomP.xtc',
		'whichframes':None,
		},
	'v532-20000-58000-100':{
		'sysname':'membrane-v532',
		'sysname_lookup':'membrane-v532-atomP',
		'director':director_aamd_asymmetric,'selector':selector_aamd_asymmetric,'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-trestles-md.part0007.20000-58000-100.atomP.xtc',
		'whichframes':None,
		},
		
	'v530-40000-90000-50':{
		'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'u5-sim-trestles-md.part0007.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Na$^{+}$',
		'ion_name':'Na',
		'composition_name':'asymmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},
	'v531-40000-90000-50':{
		'sysname':'membrane-v531',
		'sysname_lookup':'membrane-v531-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s5-sim-kraken-md.part0016.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Mg$^{2+}$',
		'ion_name':'Mg',
		'composition_name':'asymmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},		
	'v532-40000-90000-50':{
		'sysname':'membrane-v532',
		'sysname_lookup':'membrane-v532-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s5-sim-kraken-md.part0015.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Ca$^{2+}$',
		'ion_name':'Ca',
		'composition_name':'asymmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},

	'v509-40000-90000-50':{
		'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s6-kraken-md.part0018.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Na$^{+}$',
		'ion_name':'Na',
		'composition_name':'symmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},
	'v510-40000-90000-50':{
		'sysname':'membrane-v510',
		'sysname_lookup':'membrane-v510-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s8-kraken-md.part0021.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Mg$^{2+}$',
		'ion_name':'Mg',
		'composition_name':'symmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},		
	'v511-40000-90000-50':{
		'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s8-kraken-md.part0021.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'PI2P',
		'ptdins_label':'PtdIns(4,5)P$_2$',
		'ion_label':'Ca$^{2+}$',
		'ion_name':'Ca',
		'composition_name':'symmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		},

	'v533-40000-90000-50':{
		'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-kraken-md.part0015.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'P35P',
		'ptdins_label':'PtdIns(3,5)P$_2$',
		'ion_label':'Mg$^{2+}$',
		'ion_name':'Mg',
		'composition_name':'asymmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		},		
	'v534-40000-90000-50':{
		'sysname':'membrane-v534',
		'sysname_lookup':'membrane-v534-keys',
		'director':director_aamd_asymmetric,
		'selector':selector_aamd_asymmetric,
		'protein_select':None,
		'residues':'infer',
		'trajsel':'s4-sim-kraken-md.part0013.40000-90000-50.keylipidatoms.xtc',
		'whichframes':None,
		'ptdins_resname':'P35P',
		'ptdins_label':'PtdIns(3,5)P$_2$',
		'ion_label':'Ca$^{2+}$',
		'ion_name':'Ca',
		'composition_name':'asymmetric',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		},

	}
	
#---choose default midplane resolution
gridspacing = 0.5
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'
	
#---coallate two dictionaries
master_dict = dict()
for key in analysis_descriptors.keys():
	if key in analysis_descriptors_extra.keys():
		master_dict.update({key:dict(analysis_descriptors[key],
			**analysis_descriptors_extra[key])})
	else: master_dict.update({key:analysis_descriptors[key]})
analysis_descriptors = master_dict

#---COLOR HEADER
#-------------------------------------------------------------------------------------------------------------

def color_dictionary_aamd(ionname=None,lipid_resname=None,comparison=None):
	'''Master listing of colors for different comparisons, AAMD version.'''

	#---name the colors from brewer here, used as the default set for now
	clrs = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors
	clrs_names = {
		'red':0,
		'blue':1,
		'green':2,
		'purple':3,
		'orange':4,
		'yellow':5,
		'brown':6,
		'lavender':7,
		'grey':8,
		}
	
	#---compare different cation types
	if comparison == 'ions':
		colordict = {
			'Na':'green',
			'Mg':'red',
			'Ca':'blue',
			}
		return clrs[clrs_names[colordict[ionname]]]

	#---compare cations and lipid types
	if comparison == 'ions_phospate_position':
		colordict = {
			('Mg','PI2P'):'red',
			('Ca','PI2P'):'blue',
			('Mg','P35P'):'green',
			('Ca','P35P'):'purple',
			}
		return clrs[clrs_names[colordict[(ionname,lipid_resname)]]]

	#---return black
	else: return 'k'
	
def fill_between(x,y1,y2=0,ax=None,new_alpha=0.75,**kwargs):
	'''Wrap ax.fill_between so that the legend has a rectangle.'''
	ax = ax if ax is not None else plt.gca()
	ax.fill_between(x, y1, y2, **kwargs)
	fixargs = kwargs
	fixargs['alpha'] = new_alpha
	p = plt.Rectangle((0, 0), 0, 0, **fixargs)
	ax.add_patch(p)
	return p

#!/usr/bin/python

#---settings
basedir = '~/compbio'
paramfile = 'parameters.in'

#---parameter file specifications
sqlvartypes = {'float':'real','int':'int','str':'text'}
param_names = [
	'kappa','blen','gsize','sys_geom','depth','period','n_clathrin','J_Clath','c_range','c_curv','n_annulus',
	'Kzero','neps','C_0','R_2','aratio','ran_seed','memb_restart','init_file','Process_torun'
	]
param_types = {
	'kappa':'float',
	'blen':'float',
	'gsize':'int',
	'sys_geom':'str',
	'depth':'float',
	'period':'float',
	'n_clathrin':'int',
	'J_Clath':'float',
	'c_range':'float',
	'c_curv':'float',
	'n_annulus':'int',
	'Kzero':'int',
	'neps':'int',
	'C_0':'float',
	'R_2':'float',
	'aratio':'float',
	'ran_seed':'float',
	'memb_restart':'float',
	'init_file':'str',
	'Process_torun':'str',
	}
	
#---generic type converter
sqltypes = {
	str:'text',
	int:'int',
	float:'real',
	}

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------
	
def show_query(title, qry):
	'''Shortcut to display a query to standard output.'''
	status('status: QUERY: %s' % (title))
	cur.execute(qry)
	for row in cur.fetchall():
		status(row)
	status('status: done query')

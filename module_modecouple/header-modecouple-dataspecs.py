#!/usr/bin/python

import re

'''
Requires: 
	storename : the table name and the pkl prefix name
	pkl_metadata : dictionary of relevant parameters and associated types for each pickle specification
'''

storename = 'modecouple'
pkl_metadata = {
	'callsign':'text',
	'term':'text',
	'pklname':'text',
	'kappa':'float',
	'gamma':'float',
	'curvature':'text'
	}

#---? required to specify how to refresh
refresh_requirements = {
	'term':['hqs_hqps'],
	}

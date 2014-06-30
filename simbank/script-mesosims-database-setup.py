#!/usr/bin/python -i

import sys
import os
import re
import psycopg2

#---allow imports from parent directory
sys.path.insert(0,os.path.abspath('..'))

from membrainrunner import *

#---load dictionaries
execfile('header-mesosims-database.py')

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def show_query(title, qry):
	'''Shortcut to display a query to standard output.'''
	status('status: QUERY: %s' % (title))
	cur.execute(qry)
	for row in cur.fetchall():
		status(row)
	status('status: done query')

def regen_mesosims(cur):
	"""Regenerate the blank mesosims (and related) databases from parameter file specifications."""

	#---if mesosim table already exists, drop it (beware)
	cur.execute("select * from information_schema.tables where table_name=%s", ('mesosims',))
	if bool(cur.rowcount):
		status('status: mesosims table exists')
		msg = 'status: okay to drop the table?'
		shall = True if raw_input("%s (y/N) " % 'status: okay to drop the table?').lower() == 'y' else False
		confirmed = True if raw_input("%s (y/N) " % 'status: confirmed?').lower() == 'y' else False
		if shall and confirmed:
			cur.execute('DROP TABLE mesosims')
			conn.commit()

	#---if mesosim_datastore table already exists, drop it (beware)
	cur.execute("select * from information_schema.tables where table_name=%s", ('mesosims_datastore',))
	if bool(cur.rowcount):
		status('status: mesosims table exists')
		msg = 'status: okay to drop the table?'
		shall = True if raw_input("%s (y/N) " % 'status: okay to drop the table?').lower() == 'y' else False
		confirmed = True if raw_input("%s (y/N) " % 'status: confirmed?').lower() == 'y' else False
		if shall and confirmed:
			cur.execute('DROP TABLE mesosims_datastore')
			conn.commit()

	#---generate the mesosims table with a key, callsign, and a root path to that group of simulations
	cur.execute("CREATE TABLE mesosims_datastore (id serial PRIMARY KEY, callsign text, path text)")
	conn.commit()

	#---generate the mesosims table with a key, callsign, and a root path to that group of simulations
	cur.execute("CREATE TABLE mesosims (id serial PRIMARY KEY, callsign text, path text)")
	conn.commit()

	#---load the mesosims table with the appropriate parameters
	for param in param_names:
		if 0: cur.execute("ALTER TABLE mesosims ADD COLUMN %s %s" % (param,sqlvartypes[param_types[param]],))
		cur.execute('ALTER TABLE mesosims ADD COLUMN '+param+' '+sqlvartypes[param_types[param]])
		conn.commit()
	
	#---list current set of columns in mesosims table
	cur.execute("SELECT * FROM mesosims LIMIT 0")
	colnames = [desc[0] for desc in cur.description]
	print colnames
	
def pathfinder(basedir,include_patterns,file_in_path=None,excludes=None,target='filename',
	immediate=False):
	'''Search a directory tree for a particular kind of directory.'''
	#---walk the directory tree and save paths containing the parameter file
	if file_in_path != None:
		paths_with_parameters = []
		pattern = re.compile(file_in_path)
		rootdir = os.path.expanduser(basedir)
		for root,dirnames,filenames in os.walk(rootdir,followlinks=True):
			for filename in filenames:
				path = root+filename
				m = pattern.match(filename)
				if m and target == 'filename': paths_with_parameters.append(path)
				elif m and target == 'pathname': paths_with_parameters.append(root)
		paths_with_parameters = list(set(paths_with_parameters))
	else:
		if immediate:
			paths_with_parameters = os.listdir(basedir)			
		else:
			paths_with_parameters = []
			rootdir = os.path.expanduser(basedir)
			for root,dirnames,filenames in os.walk(rootdir,followlinks=True):
				for filename in filenames:
					paths_with_parameters.append(root+filename)
	#---filter paths to find valid ones
	validpaths = []
	#---excluded patterns
	exclude_pattern_list = []
	if excludes != None:
		for pat in (excludes if type(excludes) == list else [excludes]):
			exclude_pattern_list.append(re.compile(pat))
	#---included patterns
	include_pattern_list = []
	for pat in (include_patterns if type(include_patterns) == list else [include_patterns]):
		include_pattern_list.append(re.compile(pat))
	#---filter all paths
	for path in paths_with_parameters:
		if all([i.match(path) == None for i in exclude_pattern_list]) and \
			all([i.match(path) != None for i in include_pattern_list]):
			validpaths.append(path)
	return validpaths
	
def read_meso_parameters(loc):
	'''Load mesoscale simulation parameters into a dictionary.'''
	thisdict = {}
	fp = open(loc,'r')
	counter = 0
	for line in fp:
		if counter < len(param_names): 
			#---convert to float if possible
			try:
				float(line.split()[0])
				paramval = float(line.split()[0])
			except ValueError:
				paramval = line.split()[0]
			thisdict[param_names[counter]] = paramval
		counter += 1
	fp.close()
	return thisdict

#---PREPARE DATABASE
#-------------------------------------------------------------------------------------------------------------

#---connect
if 'conn' not in globals(): 
	try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	except: print "I am unable to connect to the database"
	cur = conn.cursor()
	#---show a list of databases
	show_query('available databases', 'SELECT * FROM pg_database')

	#---regenerate the mesosims database
	regen_mesosims(cur)

#---PARSE DIRECTORY
#-------------------------------------------------------------------------------------------------------------

if 'validpaths' not in globals(): 
	validpaths = pathfinder(basedir,'.+/mesomembrane-v[0-9]{4}/.+',
		excludes=['.+serial-code','.*/run_settings\.[0-9]{4}'],file_in_path=paramfile,target='pathname')
	pklnames = pathfinder('/home/rpb/worker/repo-pickles/','.+meso.+\.pkl',immediate=True)

#---LOAD INTO DATABASE
#-------------------------------------------------------------------------------------------------------------

#---select a single entry
for simpath in validpaths:
	status('status: parsing '+simpath)

	#---compile parameters and path
	thisdict = read_meso_parameters(simpath+'/'+paramfile)
	thisdict['path'] = simpath
	thisdict['callsign'] = simpath.split('mesomembrane-')[1].split('/')[0]

	#---insert data
	cur.execute('INSERT INTO mesosims ('+\
		','.join(['callsign','path']+param_names)+') VALUES ('+\
		', '.join(["%("+i+")s" for i in ['callsign','path']+param_names])+');',
		thisdict)
	conn.commit()


'INSERT INTO mesosims ('+','.join(['callsign','path']+param_names)+') VALUES ('+', '.join(["%("+i+")s" for i in ['callsign','path']+param_names])+');'


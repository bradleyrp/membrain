#!/usr/bin/python

'''
Construct the tables necessary for the bookkeeping for one calculation.\n
Note that this is a one-time use database.

PSUEDOCODE
	1. get column names incoming from the script that instantiates this class
	2. make a table corresponding to pickle objects
	3. If this is a three-way calculation
		import functions for parsing the simulations
		parse the simulations
		use the results to load up other DOWNSTREAM tables
			i.e. the table for storing kappa values etc "mesosims_datastore" after the structure is 
				...computed and the simulation is identified
		
'''

import sys,os,re
import psycopg2
import psycopg2.extras

sys.path.insert(0,os.path.abspath('..'))
from membrainrunner import *

#---CLASS
#-------------------------------------------------------------------------------------------------------------

class DataFace:

	def __init__(self,**kwargs):
		'''
		Initiate a connection to the SQL database.
		'''
		
		#---catch settings
		bigkeylist = ['dbconnect_string','pickles']
		self.sets = dict()
		for key in bigkeylist: self.sets[key] = kwargs[key]
		
		#---connect
		try: self.conn = psycopg2.connect(self.sets['dbconnect_string'])
		except: raise Exception('except: unable to connect to database')
		self.cur = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
	
	def ec(self,cmd):
		'''Execute and commit a change to the database.'''
		self.cur.execute(cmd)
		self.conn.commit()

	def query(self,qry):
		'''Shortcut to return a database query.'''
		self.cur.execute(qry)
		return self.cur.fetchall()

	def scan_pklfiles(self):
		'''Return a list of files in the datastore (repo-pickles).'''
		rootdir = os.path.expanduser(self.sets['pickles'])
		pklfilenames = []
		for (dirpath, dirnames, filenames) in os.walk(rootdir):
			pklfilenames += filenames
			break
		return pklfilenames

	def unique(self,specs,extras=None):
		'''
		Function for looking up uniquely-named data files from a dictionary.
		'''
		
		#---mimics the new function which reduces a dictionary-within-a-dictionary to a string
		if type(extras) == dict:
			compact = dict()
			#---iterate one level to see if the extraspecs dictionary contains any dictionaries
			for key in extras:
				if type(extras[key]) == dict:
					compact[key] = '|'.join([subkey+':'+str(extras[key][subkey]) 
						for subkey in extras[key].keys()])
				else: compact[key] = str(extras[key])
			rowspecs = dict(compact.items()+specs.items())
		else: rowspecs = dict(specs)
		if 'i' in rowspecs.keys(): 
			rowspecs['id'] = rowspecs['i']
			del rowspecs['i']
		if 'calculation' in rowspecs.keys(): del rowspecs['calculation']
		#---check the database
		lookup = self.query('SELECT * FROM dataref_modecouple '+\
			' WHERE '+' AND '.join(['( '+key+'=\''+str(rowspecs[key])+'\' OR '+key+' IS NULL)' 
			for key in rowspecs.keys()]))
		if lookup == []: return None
		elif len(lookup) > 1: raise Exception('multiple entries match this query')
		else: return lookup[0]
		
	def pklsig(self,pklname):

		'''
		Function which decomposes a pickle file name into metadata required for entry into the database.\n
		Note that the filename may not code all of the metadata required to understand the object.
		In this case, the filename metadata must have a unique key necessary to pull the rest of the 
		information from the database. This should be assured by the namer function.
		'''

		sigs = dict()
		#---get callsign
		callsign = pklname.split('.')[2]
		if not re.match('^v[0-9]{3,4}\-[0-9]+\-[0-9]+\-[0-9]+$',callsign):
			raise Exception('except: invalid callsign in pickle file name')
		else: sigs['callsign'] = callsign
		#---set pickle file name
		sigs['pklname'] = pklname
		#---get other signifiers, which must come in key,value pairs separated by hyphens
		#---...located within the fourth dot-delimited section of the filename
		siglist = pklname.split('.')[3].split('-')
		if len(siglist)%2 != 0: 
			raise Exception('except: invalid number of extra signifiers in the filename')
		else:
			for i in range(len(siglist)/2): 
				sigs[siglist[i*2]] = siglist[i*2+1]
		return sigs
		
	def namer(self,specs,index=None):
		'''Makes a pickle name from specs and possibly an index if the database records other metadata.'''
		#---data type (calculation) and callsign are listed first by convention
		basename = 'pkl.'+specs['calculation']+'.'+specs['callsign']
		for key in ['calculation','callsign']:
			if key in specs.keys(): del specs[key]
		basename += ('.' if len(specs.keys())>0 else '')
		basename += '-'.join([key+'-'+specs[key] for key in specs.keys()])
		if type(index) == int: basename += '-i-'+str(index)
		return basename
		
	def confirm_refresh(self,sigs):
		'''
		Confirm that this particular pickle can be added to the database.\n
		Note that this is function is for development only.
		'''
		if 'refresh_requirements' in self.sets.keys():
			for key in self.sets['refresh_requirements'].keys():
				if sigs[key] not in self.sets['refresh_requirements'][key]: return False
		return True
				
	def refresh_dataref(self,**dataspecs):

		'''
		Function which updates the dataref table for a particular calculation 
		by scanning for valid pickles.\n
		Note that this function handles a large part of the pickle-to-database interface.
		'''

		#---catch important settings
		self.storename = dataspecs['storename']
		bigkeylist = ['pkl_metadata','refresh_requirements']
		for key in bigkeylist: 
			if key in dataspecs.keys():
				self.sets[key] = dataspecs[key]
		
		#---make the table if necessary
		self.cur.execute("select * from information_schema.tables where table_name=%s",
			('dataref_'+self.storename,))
		if not bool(self.cur.rowcount): 
			status('status: building dataref_'+self.storename+' table')
			self.make_dataref() 

		#---find relevant pickles according to their prefix which must follow "pkl." in the filename
		pkllist = [fn for fn in self.scan_pklfiles() if re.match('^pkl.'+self.storename+'.+',fn)]
		#---extract further signifiers from the filename according to simple rules
		#---...that is, the second dot delimited feature should be the call sign
		#---...then all subsequent parameters are explicitly listed as hyphen-delimited key-data pairs
		for pklname in pkllist:
			#---get metadata from the filename
			sigs = self.pklsig(pklname)
			#---add entry to the database if necessary
			if 'i' in sigs.keys(): 
				sigs['id'] = sigs['i']
				del sigs['i']
			self.ec('SELECT * FROM dataref_'+self.storename+' WHERE '+\
				' AND '.join([key+'=\''+sigs[key]+'\'' for key in sigs.keys()]))
			if not bool(self.cur.rowcount) and self.confirm_refresh(sigs):
				self.cur.execute('INSERT INTO dataref_'+self.storename+' ('+\
					','.join(sigs.keys())+') VALUES ('+\
					', '.join(["%("+i+")s" for i in sigs.keys()])+');',
					sigs)
				self.conn.commit()
			del sigs
			
		#---perform the check in the opposite direction to confirm database pkl files still exist
		#---remove entries for missing pickle files
		self.cur.execute('SELECT * FROM dataref_'+self.storename)
		rows_to_delete = []
		for row in self.cur: 
			if row['pklname'] == None or not os.path.isfile(self.sets['pickles']+row['pklname']):
				status('error: missing pklname/id = '+\
					(str(row['id']) if row['pklname'] == None else row['pklname']))
				status('status: removing row ---> '+str(row))
				rows_to_delete.append(row['id'])
		if rows_to_delete != []: raw_input('warning: anything to continue but beware row deletions')
		for rowid in rows_to_delete:
			self.ec('DELETE FROM dataref_'+self.storename+' WHERE id='+str(rowid))
			
	def new(self,specs,extras=None):
	
		'''
		Core database functionality in which generates a new entry with excess metadata.
		'''

		if type(extras) == dict:
			compact = dict()
			#---iterate one level to see if the extraspecs dictionary contains any dictionaries
			for key in extras:
				if type(extras[key]) == dict:
					compact[key] = '|'.join([subkey+':'+str(extras[key][subkey]) 
						for subkey in extras[key].keys()])
				else: compact[key] = str(extras[key])
			newrow = dict(compact.items()+specs.items())
		else: newrow = dict(specs)
		if 'calculation' in newrow.keys(): del newrow['calculation']
		self.cur.execute('INSERT INTO dataref_'+self.storename+' ('+\
			','.join(newrow.keys())+') VALUES ('+\
			','.join(["%("+i+")s" for i in newrow.keys()])+') RETURNING id;',
			newrow)
		index = self.cur.fetchone()[0]
		self.conn.commit()
		return index
		
	def update(self,ind,**kwargs):
		'''Update an entry.\n
		This function is designed primarily to insert the name of a data file that was created for a 
		pre-existing record.'''
		self.ec('UPDATE '+kwargs['table']+' SET ('+\
			','.join([key for key in kwargs.keys() if key != 'table'])+') = ('+\
			','.join(['\''+kwargs[key]+'\'' for key in kwargs.keys() if key != 'table'])+\
			') WHERE id=\''+str(ind)+'\';')

	def make_dataref(self):

		'''
		Generate a "dataref" table in the database to store connections to pickle files in the repository.\n
		This function uses the storename passed to refresh_dataref to name the table while the column names 
		and corresponding types are provided via pkl_metadata defined in the calculation-specific header file
		which is also passed via keyword arguments to refresh_dataref and then filtered and stored as 
		self.sets for use here.
		'''

		#---create the table with columns
		status('status: creating table')
		cmd = 'CREATE TABLE '+('dataref_'+self.storename)+' (id serial PRIMARY KEY'
		for key in (self.sets['pkl_metadata']): cmd += ', '+key+' '+(self.sets['pkl_metadata'])[key]
		self.ec(cmd+')')
		
		

#!/usr/bin/python

import sys,os,re,time
import datetime
import psycopg2
import psycopg2.extras

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---settings
lenscale = 0.5
startframe = 2000
endframe = 5000
calctypes = ['structure','update_kappas'][-1:]

#---INCLUDES
#-------------------------------------------------------------------------------------------------------------

#---allow imports from parent directory
sys.path.insert(0,os.path.abspath('..'))

#---initiate lab notebook
loggerdir = '/home/rpb/worker/repo-journals'+'/'
logname = 'labnotes-calculation-'+'-'.join(calctypes)+'-time-'+\
	datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M.%S')
if not os.path.exists(loggerdir+logname): os.makedirs(loggerdir+logname)
else: raise Exception('except: journal directory exists')
cwd = os.path.realpath('.')+'/'
files_to_log = ['header-mesosims-database.py',sys.argv[0].strip('./')]
for f in files_to_log: os.system('cp '+cwd+f+' '+loggerdir+'/'+logname+'/')
logfile = loggerdir+logname+'/log-'+logname

from membrainrunner import *
execfile('../locations.py')

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---load dictionaries
execfile('header-mesosims-database.py')
	
#---list associated pickle files
assoc_pickles = {
	'structure':{
		'callsign':'text',
		'path':'text',
		'parent_mesosims':'int',
		'rep':'int',
		'startframe':'int',
		'endframe':'int',
		'lenscale':'real',
		'simtype':'text',
		'locate':'text',
		'nbase':'int',
		'shortname':'text',
		'hascurv':'boolean',
		'testname':'text',
		'rundirnum':'int',
		}
	}

#---suffixes for files in the datastore (repo-pickles)
filename_suffixes = ['.png','.pkl']
name_params = ['C_0','rundir']

#---lookup_logic gives parameters which must match when comparing dictionaries of metadat
#---list of lists, each of which must have at least one matching element between dicts
lookup_logic = [['C_0','rundir'],['callsign']]
lookup_logic = [['rundir'],['callsign']]

#---UTILITY FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def filename_parser(filename,substring,prefix):
	'''Checks a filename for a numerical parameter and returns it.'''
	if re.search(substring,filename) and re.search(prefix+'.',filename):
		chopstring = filename
		for suf in filename_suffixes:
			chopstring = chopstring.strip(suf)
		fnlist = chopstring.strip('pkl.').strip('pkl').split('-')
		ind = fnlist.index(substring)
		return float(fnlist[ind+1])
	else: return None

def build_dict(cursor,row):
	'''Load a row into a dictionary from postgresql.'''
	thisdict = {}
	for key,col in enumerate(cursor.description): thisdict[col[0]] = row[key]
	return thisdict
	
def scan_pklfiles():
	'''Return a list of files in the datastore (repo-pickles).'''
	rootdir = os.path.expanduser(pickles)
	pklfilenames = []
	for (dirpath, dirnames, filenames) in os.walk(rootdir):
		pklfilenames += filenames
		break
	return pklfilenames
	
def pickle_categorizer(pklname,dattype):
	'''
	Takes a pickle file name and makes intelligent guesses about what the file contains.\n
	Note that this procedure handles a fairly fluid but still specific set of naming conventions.\n
	'''
	metadat = dict()
	#---filter by pickle type
	#---? needs generalized
	if pklname.split('.')[1] == 'structures' and pklname.split('.')[2] == 'meso':
		#---get callsign
		callsign = [i for i in pklname.replace('-','.').split('.') if i[0] == 'v']
		callsign = callsign[0] if len(callsign) > 0 else None
		metadat['callsign'] = callsign
		for name in name_params:
			val = filename_parser(pklname,name,'pkl')
			if val != None: metadat[name] = val		
		return [pklname,metadat]
	else: return None

def collect_extra_metadata(rowdict):
	'''This function collects additional metadata about a particular row in mesosims.'''
	#---check for rundir, assumes directory follows format "rundir-%d"
	if 'path' in rowdict.keys():
		path = rowdict['path']
		if re.search('rundir',path):
			modifier = path.split('rundir-')[1]
			rowdict['rundir'] = int(modifier)
	#---check for replicate number, assumes directory follows format "rep-%d"
	if 'path' in rowdict.keys():
		path = rowdict['path']
		if re.search('rep',path): 
			modifier = path.split('rep-')[1]
			rowdict['rep'] = int(modifier)
	return rowdict

def scan_datastore(dataref_key):
	'''
	The scan function parses the datastore, specifically pickle files in repo-pickles, for elements that
	correspond to entries in the simulations database table called "mesosims". The results are added to the
	associated table.
	'''
	#---inner function compares metadata and row dictionaries to see if a pickle matches a database row
	#---excessive use of lowercase function is due to the fact that postgresql is case insensitive
	inner = lambda r,m,i: (
		(i in r.keys() or i.lower() in r.keys()) and 
		(i in m.keys() or i.lower() in m.keys()) and 
		(r[(i if i in r.keys() else i.lower())] == m[(i if i in m.keys() else i.lower())]))
	#---batch categorizer loops over all pickle files and extracts a dictionary of filename attributes
	valids = []
	for pklfile in scan_pklfiles(): 
		cat = pickle_categorizer(pklfile,'structmeso') #---? needs generalized
		if cat != None: valids.append(cat)
	#---intable holds pickles that found a match in the mesosims table
	intable = []
	#---outtable holds pickles that could not be matched with a row in the mesosims table
	outtable = []
	#---loop over pickles to find associated simulations
	for pklname,metadat in valids:
		metadat['pklfile'] = pklname
		#---filter by callsign
		cur.execute('SELECT * FROM mesosims')
		match = False
		for row in cur:
			rowdict = dict(row)
			rowdict = collect_extra_metadata(rowdict) #---? better naming
			#---parse the lookup_logic list to see if the dictionaries match
			if all([any([inner(metadat,rowdict,l2) for l2 in l]) for l in lookup_logic]):
				intable.append([metadat,rowdict])
				match = True
				break
		if not match: outtable.append(metadat) #---? still useful to add outtable?
	return intable,outtable

#---DATAREF DATABASE FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def create_dataref_table(dataref_key):	
	'''Adds a new pointer to data in the datastore (repo-pickles) for a particular type of data.'''
	tablename = 'dataref_'+dataref_key
	#---check if table exists
	cur.execute("select * from information_schema.tables where table_name=%s", (tablename,))
	#---delete table if it exists (TEMPORARY, however this is very fast)
	if bool(cur.rowcount): 
		status('status: '+tablename+' table exists')
		#msg = 'status: okay to drop the table?'
		#shall = True if raw_input("%s (y/N) " % 'status: okay to drop the table?').lower() == 'y' else False
		#confirmed = True if raw_input("%s (y/N) " % 'status: confirmed?').lower() == 'y' else False
		#---override
		shall,confirmed = True,True
		if shall and confirmed:
			cur.execute('DROP TABLE '+tablename)
			conn.commit()
	#---generate the table with columns specified by global assoc_pickles dictionary
	cmd = 'CREATE TABLE '+('dataref_'+dataref_key)+' (id serial PRIMARY KEY'
	for key in (assoc_pickles[dataref_key]):
		cmd += ', '+key+' '+(assoc_pickles[dataref_key])[key]
	cmd += ')'
	cur.execute(cmd)
	conn.commit()
	
def load_dataref_table(dataref_key,intable):
	'''Load a list of pairs of metadata for pickles and rows in mesosims into a dataref table.'''
	tablename = 'dataref_'+dataref_key
	#---compute a set of necessary columns
	#---note that assoc_pickles has some defaults already
	newcolumns = []
	for metadat,rowdict in intable:
		for i in metadat.keys():
			if i not in [k[0] for k in newcolumns]:
				newcolumns.append([i,sqltypes[type(metadat[i])]])
	cur.execute('SELECT * from '+tablename)
	colnames = [desc[0] for desc in cur.description]
	#---add columns to the (presumably blank) dataref table
	for col in newcolumns:
		if col[0] not in colnames:
			cur.execute('ALTER TABLE '+tablename+' ADD COLUMN '+col[0]+' '+col[1])
			conn.commit()
	for metadat,rowdict in intable:
		metadat['parent_mesosims'] = rowdict['id']
		#---extra metadata added here
		metadat['startframe'] = startframe
		metadat['endframe'] = endframe
		metadat['lenscale'] = lenscale
		metadat['simtype'] = 'meso'
		cur.execute('SELECT path from mesosims where mesosims.id='+str(rowdict['id']))
		path = cur.fetchall()[0][0]
		metadat['locate'] = path
		cur.execute('SELECT gsize from mesosims where mesosims.id='+str(rowdict['id']))
		metadat['nbase'] = cur.fetchall()[0][0]
		#---invalid because uppercase
		if 0 and 'c_0' in metadat.keys(): metadat['C_0'] = metadat['c_0']
		metadat['shortname'] = '$\mathrm{meso,}\:C_0='+'{0:.4f}'.format(metadat['C_0'])+'a_0^{-1}$'
		metadat['hascurv'] = True if metadat['C_0'] > 0 else False
		metadat['testname'] = metadat['callsign']+'-C_0-'+str(metadat['C_0'])
		#---if no explicit replicate number, add a null value
		if 'rep' not in metadat.keys(): metadat['rep'] = -1	
		#---insert data
		cur.execute('INSERT INTO '+tablename+' ('+\
			','.join(metadat.keys())+') VALUES ('+\
			', '.join(["%("+i+")s" for i in metadat.keys()])+');',
			metadat)
		conn.commit()

def missing_calcs(dataref):
	'''Find rows in mesosims with no corresponding entry in a dataref table.'''
	tablename = 'dataref_'+dataref
	cur.execute('SELECT * FROM mesosims WHERE id NOT IN (SELECT mesosims.id FROM mesosims,'+\
		tablename+' WHERE mesosims.id='+tablename+'.parent_mesosims)')
	return [dict(i) for i in cur.fetchall()]
	
def calculate_structures(paramdict,start,end,lenscale):
	'''Calculate structure pickles for mesoscale simulations.'''
	#---loop over needed pickles
	for newcalc in paramdict:
		#---check for replicate subfolders or continue if using rundir convention
		rootdir = newcalc['path']
		repdirs = []
		for (dirpath, dirnames, filenames) in os.walk(rootdir):
			repdirs += dirnames
			break
		params = dict(newcalc)
		if any([i[:4]=='rep-' for i in dirnames]): 
			raise Exception('except: not prepared for replicate subdirectories')
		#---the following section disables this code for use on directory structures not organized by rundirs
		if params['path'].split('/')[-1].split('-')[0] != 'rundir':
			raise Exception('except: only set for rundir-style directories')
		#---augment the dictionary with useful information
		#---note some of these are redundant or less general but this keeps with previous conventions
		#---the following parameters were moved to the database via load_dataref_table and disabled below
		if 0:
			#---note: changed "end" to "endframe" because it's a keyword in sql
			params['startframe'] = startframe
			params['endframe'] = endframe
			params['lenscale'] = lenscale
			params['simtype'] = 'meso'
			params['locate'] = params['path']
			params['nbase'] = params['gsize']
			params['shortname'] = r'$\mathrm{meso,}\:C_0=0.014a_0^{-1}$'
			if 'c_0' in params.keys(): params['C_0'] = params['c_0']
			params['hascurv'] = True if params['C_0'] > 0 else False
			params['testname'] = params['callsign']+'-C_0-'+str(params['C_0'])
		#---enabled the lenscale when testing the code on older pickles not stored in rundir
		if 0: params['lenscale'] = lenscale
		params['timestamp_created'] = \
			datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M.%S')
		rundirnum = int(params['path'].split('/')[-1].split('-')[1])
		pklname = 'pkl.structures.meso.'+\
			params['callsign']+'-'+\
			'C_0-'+str(params['C_0'] if 'C_0' in params.keys() else params['c_0'])+'-'+\
			'len-'+str(params['lenscale'])+'-'+\
			'rundir-'+str(rundirnum)+\
			'.pkl'
		mset = unpickle(pickles+pklname)
		status('status: checking for structure pickle, '+params['testname'])
		if mset == None:
			status('status: calculation = '+pklname)
			mset = MembraneSet()
			status('status: computing structure, '+params['testname'])
			c0sraw = array(mset.load_points_vtu(params['locate'],extra_props='induced_cur',
				start=params['startframe'],end=params['endframe'],nbase=params['nbase'],
				lenscale=lenscale,prefix='EQUIB-conf-'))[:,0]
			mset.surfacer()
			c0s = mset.surfacer_general(c0sraw)
			result_data = MembraneData('c0map')
			result_data.data = array(array(c0s))
			for i in params: 
				result_data.addnote([i,params[i]])
			mset.store.append(result_data)
			pickledump(mset,pklname,directory=pickles)
			
def update_kappas():
	'''Add observed bending rigidity values to the database.'''	
	#---check for kappa_apparent column in dataref_structure
	tablename = 'dataref_structure'
	cur.execute('select * from '+tablename)
	colnames = [desc[0] for desc in cur.description]
	if 'kappa_apparent' not in colnames:
		cur.execute('ALTER TABLE '+tablename+' ADD COLUMN kappa_apparent real')
		conn.commit()
	#---get pickle files for unknown kappa_apparent values
	cur.execute('SELECT pklfile FROM dataref_structure WHERE kappa_apparent IS NULL')
	pkllist = [i[0] for i in cur.fetchall()]
	#---loop over needed pickles
	for pklfile in pkllist:
		mset = unpickle(pickles+pklfile)
		if mset.undulate_kappa == 0.0: mset.calculate_undulations()
		status('status: checking for structure pickle '+pklfile)
		print 'kappa = '+str(mset.undulate_kappa)
		cur.execute("select id from dataref_structure where pklfile='"+pklfile+"'")
		rowid = cur.fetchall()[0][0]
		cur.execute('UPDATE '+tablename+' SET kappa_apparent=%s where id=%s',(mset.undulate_kappa,rowid,))
		conn.commit()

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---connect
if 'conn' not in globals(): 
	try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	except: raise Exception('except: cannot connect to database')
	cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

#---make the table, scan for structure pickles, populate the dataref_structure table
for calctype in calctypes:
	if not re.search('update_',calctype):
		create_dataref_table(calctype)
		intable,outtable = scan_datastore(calctype)
		load_dataref_table(calctype,intable)
		needs = missing_calcs(calctype)

#---settings
if 'structure' in calctypes: needs = [i for i in needs if int(i['callsign'][1:]) == 2014]

#---computations
if 'structure' in calctypes: calculate_structures(needs,startframe,endframe,lenscale)
if 'update_kappas' in calctypes: update_kappas()

"""
NOTE: extremely useful way to lookup items in a child table and print alongside items from the parent:

SELECT kappa_apparent,parent_mesosims,mesosims.kappa 
FROM dataref_structure,mesosims 
WHERE kappa_apparent IS NOT NULL and mesosims.id=dataref_structure.parent_mesosims;
"""



	
	

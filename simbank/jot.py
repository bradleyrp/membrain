#!/usr/bin/python

dataref_key = 'structure'
if 1:
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
tmp = [i for i in valids if i[1]['callsign']=='v2014']

#!/usr/bin/python

#---connect
if 'conn' not in globals(): 
	conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	#try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	#except: raise Exception('except: cannot connect to database')
	cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
	
#---retrieve a set of experiments for comparison
cur.execute('SELECT * from mesosims')
select_mesosims = [dict(i) for i in cur.fetchall()]
cur.execute('SELECT * from dataref_structure')
select = [dict(i) for i in cur.fetchall()]
#---note that I am filtering the table in python and not postgresql
select = [i for i in select if all([
	i[j] == select_criteria_meso[j] for j in select_criteria_meso.keys()])]
#---populate analysis descriptors from the database
for params in select:
	ind = where([i['id']==params['parent_mesosims'] for i in select_mesosims])[0][0]
	combo = dict(params.items() + select_mesosims[ind].items())
	rundirnum = int(combo['rundir'])
	key = params['callsign']+'-rundir-'+str(rundirnum)
	#---always fix uppercase naming when importing to python
	if 'c_0' in combo.keys(): combo['C_0'] = combo['c_0']
	combo['detail_name'] = combo['shortname']
	analysis_descriptors[key] = combo

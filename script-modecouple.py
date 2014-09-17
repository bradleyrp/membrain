#!/usr/bin/python -i

from membrainrunner import *
allsets = dict()
execfile('locations.py',allsets)
execfile('locations/header.py',allsets)
from module_modecouple.ModeCouple import *
from module_modecouple.ModeSum import *
from module_database.DataFace import DataFace
import copy

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---select a simulation system
callsign = [
	'v550-300000-400000-200',
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	][-1]
	
#---hypothesis for the calculation section which is only a subset of the full hypothesis
hypothesis_calc = {
	'curvature':{
		'type':'dimple',
		'C_0':0.001,
		'sigma_a':10,
		'sigma_b':10,
		}
	}
	
#---the default hypothesis around which we sweep parameters
hypothesis_default = {
	'curvature':{
		'type':'dimple',
		'C_0':0.001,
		'sigma_a':10,
		'sigma_b':10,
		},
	'kappa':{
		'type':'disc',
		'back':20,
		'fore':24,
		'radius':20,
		},
	'gamma':0.0,
	}

#---sweep over curvatures
sweep_curv = {
	'curvature':{
		'C_0':[0.001,0.005,0.01,0.018,0.02,0.022,0.024,0.03,0.035,0.04,0.05]
		},
	}

#---sweep over bending rigidity
sweep_kappa = {
	'kappa':{
		'fore':[20,22,24,28,32],
		},
	}
	
#---combine sweeps
sweep = dict(sweep_curv.items()+sweep_kappa.items())
	
#---construct hypotheses from sweep variables and the default hypothesis
hypotheses = []
for topkey in sweep.keys():
	if type(sweep[topkey]) == dict:
		for key in sweep[topkey].keys():
			for val in list(sweep[topkey][key]):
				newhypo = copy.deepcopy(hypothesis_default)
				newhypo[topkey][key] = val
				hypotheses.append(newhypo)
				del newhypo
	else:
		for val in sweep[topkey]:
			newhypo = copy.deepcopy(hypothesis_default)
			newhypo[topkey] = val
			hypotheses.append(newhypo)
			del newhypo

#---choose what to do
routine = [
	'calc',
	'hypothesize',
	][1:]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def specfilter_calc(specs):
	''' '''
	filt = dict()
	for key in specs.keys():
		if key in hypothesis_calc.keys():
			filt[key] = specs[key]
	return filt
	

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'calc' in routine:

	'''
	Since parts of the hypothesis, such as the bending rigidity field, are easy to compute, the calculation
	section doesn't require a full hypothesis. Instead, it computes the terms which contribute to the final
	Helfrich sum and depend on a minimal hypothesis, in this case consisting only of the curvature field. 
	These results are then stored in hd5f binaries which are unpacked during the "hypothesize" stage.
	'''
	
	for hypothesis in hypotheses[:1]:
	
		status('status: INVESTIGATING HYPOTHESIS')
		status('status: '+str(hypothesis))

		#---create an interface to the database
		df = DataFace(**allsets)
		#---pass along the calculation-specific specs
		dataspecs = dict()
		execfile('./module_modecouple/header-modecouple-dataspecs.py',dataspecs)
		#---refresh the list of available pickles
		df.refresh_dataref(**dataspecs)

		#---find or compute the required terms and save to the datastore
		ms = ModeSum(hypothesis,callsign,**allsets)
		for termname in ['hqs_hqps','hqs_c0qps','c0qs_hqps','c0qs_c0qps']:
			specifier = {
				'term':termname,
				'callsign':callsign,
				'calculation':'modecouple',
				}
			#---only add hypothesis to database if the term depends on it
			if df.unique(specifier,extras=(None if termname == 'hqs_hqps' 
				else specfilter_calc(hypothesis))) == None:
				ind = df.new(specifier,
					extras=(None if termname == 'hqs_hqps' else specfilter_calc(hypothesis)))
				pklname = df.namer(specifier,index=ind)
				if ms.compute_modesum_term(term=termname,pklname=pklname):
					df.update(ind,pklname=pklname+'.h5pkl',table='dataref_modecouple')
			else: status('status: found the term '+termname+' in the repository')
		if 'ms' in globals(): del ms
		df.refresh_dataref(**dataspecs)

#---HYPOTHESIS TESTING
#-------------------------------------------------------------------------------------------------------------

if 'hypothesize' in routine:

	#---select a hypothesis
	hypothesis = hypotheses[0]

	#---create an interface to the database
	df = DataFace(**allsets)
	#---pass along the calculation-specific specs
	dataspecs = dict()
	execfile('./module_modecouple/header-modecouple-dataspecs.py',dataspecs)
	#---refresh the list of available pickles
	df.refresh_dataref(**dataspecs)

	#---find or compute the required terms and save to the datastore
	ms = ModeSum(hypothesis,callsign,**allsets)
	
	#---load computed terms
	termlist = []
	for termname in ['hqs_hqps','hqs_c0qps','c0qs_hqps','c0qs_c0qps']:
		specifier = {
			'term':termname,
			'callsign':callsign,
			'calculation':'modecouple',
			}
		status('status: unpacking '+termname)
		row = df.unique(specifier,extras=specfilter_calc(hypothesis))
		pklname = df.namer(specifier,index=row['id'])+'.h5pkl'
		termlist.append(unbinary(allsets['pickles']+pklname))
	
	#---compute the sum
	ms.summer(
		hqs_hqps=termlist[0],
		hqs_c0qps=termlist[1],
		c0qs_hqps=termlist[2],
		c0qs_c0qps=termlist[3])
		
	plt.imshow(real(ms.kqs).T,interpolation='nearest',origin='lower');plt.show()
		
	#---continue
	status('status: continue with script-modecouple-plot.py')













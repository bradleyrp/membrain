#!/usr/bin/python

analysis_descriptors = {
	'v614-120000-220000-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'testname':'v614',
		'locate':'pkl.structures.space10A.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl',
		'start':None,
		'end':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True
		},
	'v616-210000-310000-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'testname':'v614',
		'locate':'pkl.structures.space10A.membrane-v616.u6-trestles.md.part0005.210000-310000-200.pkl',
		'start':None,
		'end':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True
		},
	}
	
#---table of contents for mesoscale parameter sweeps
meso_expt_toc = {
	'v2004' : {
		'basedir':'/home/rpb/compbio/mesomembrane-v2004/',
		'R_2':25,
		'start':400,
		'end':900,
		},
	'v2005' : {
		'basedir':'/home/rpb/compbio/mesomembrane-v2005/',
		'R_2':25,
		'start':1500,
		'end':3500,
		},
	'v2006' : {
		'basedir':'/home/rpb/compbio/mesomembrane-v2006/',
		'R_2':25,
		'start':1500,
		'end':3500,
		},
	'v2008' : {
		'basedir':'/home/rpb/compbio/mesomembrane-v2008/',
		'R_2':25,
		'start':1500,
		'end':3500,
		},
	}

#---load the results of a parameter sweep into the dictionary
for key in meso_expt_toc.keys():
	basedir = (meso_expt_toc[key])['basedir']
	start = (meso_expt_toc[key])['start']
	end = (meso_expt_toc[key])['end']
	pname = subprocess.check_output('grep \"ID=\" '+basedir+'/mesomania.sh | tail -1',
		shell=True).split('\"')[1]
	paramstring = subprocess.check_output(\
		'grep \"param_list=\" '+basedir+'/mesomania.sh | tail -1',shell=True)
	param_precision = len(paramstring.strip().replace('\'','').split(' ')[1:-1][0].split('.')[-1])
	params = [float(i) for i in paramstring.strip().replace('\'','').split(' ')[1:-1]]
	(meso_expt_toc[key])['parameter_name'] = pname
	(meso_expt_toc[key])['parameter_sweep'] = params
	(meso_expt_toc[key])['basedir'] = basedir
	(meso_expt_toc[key])['param_precision'] = param_precision
	for param in params:
		analysis_descriptors[key+'-'+pname+'-'+str(param)] = {
			'simtype':'meso',
			'label':r'$\mathrm{meso}$',
			'detail_name':r'$\mathrm{meso,}\:C_0='+\
				(('{0:.'+str(param_precision)+'f}').format(param))+'a_0^{-1}$',
			'testname':key+'-'+pname+'-'+\
				(('{0:.'+str(param_precision)+'f}').format(param)),
			'locate':basedir+'exec/'+pname+'-'+\
				(('{0:.'+str(param_precision)+'f}').format(param))+'/rep-0/',
			'start':start,
			'end':end,
			'nbase':22,
			'hascurv':True,
			}

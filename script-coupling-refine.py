#!/usr/bin/python -i

interact = True
from membrainrunner import *
execfile('locations.py')

import sys,os,re,time
import datetime
import psycopg2
import psycopg2.extras
from scipy.optimize import leastsq

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---available functions
routine = [
	'compare_parameters_grid',
	'compare_parameters_grid_sweep_plot',
	'sweep_parameters_cgmd',
	'cgmd_test_newfield',
	][2:3]

#---selection criteria
select_criteria_meso = {'callsign':'v2016'}
cgmd_avail = [
	'v550-300000-400000-200',
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	]

#---CGMD simulations placed here for now
analysis_descriptors = {
	'v614-120000-220000-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'testname':'v614',
		'locate':'pkl.structures.space10A.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl',
		'startframe':None,
		'endframe':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True,
		'qmagfilter':[10**-10,0.4],
		},
	'v616-10000-209900-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'testname':'v616',
		'locate':'pkl.structures.sharp.space2A.membrane-v616.u6-trestles.md.part0004.110000-209900-200.pkl',
		'startframe':None,
		'endframe':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True,
		'qmagfilter':[10**-10,0.4],
		},
	'v616-210000-310000-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'testname':'v616',
		'locate_old':'pkl.structures.sharp.space2A.membrane-v616.u6-trestles.md.part0005.210000-310000-200.pkl',
		'locate':'pkl.structures.space10A.membrane-v616.u6-trestles.md.part0005.210000-310000-200.pkl',
		'startframe':None,
		'endframe':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True,
		'qmagfilter':[10**-10,0.4],
		},
	'v550-300000-400000-200': {
		'simtype':'md',
		'label':r'$\mathrm{control}$',
		'detail_name':r'$\mathrm{control}$',
		'testname':'v550',
		'locate':'pkl.structures.space10A.membrane-v550.s0-trajectory-full.md.part0006.300000-400000-200.pkl',
		'startframe':None,
		'endframe':None,
		'nbase':None,
		'hascurv':True,
		'hypo':['',1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True,
		'qmagfilter':[10**-10,0.4],
		},
	}
	
#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---note that these functions were pulled directly from script-coupling-adv-batch.py
#---note that we may wish to streamline this matching procedure from the script-coupling-adv-batch.py version

def fftwrap(dat,redundant=1):
	'''This function wraps the standard discrete FFT for a system with possible-redundant rows.'''
	trim = -1 if redundant == 1 else None 
	return fft.fftshift(fft.fft2(array(dat)[:trim,:trim]))

def ifftwrap(dat,redundant=1):
	'''This function wraps the standard discrete IFFT for a system with possible-redundant rows.'''
	trim = -1 if redundant == 1 else None 
	return fft.fftshift(fft.ifft2(array(dat)[:trim,:trim]))
	
def autocorr(dat,direct=0,lenscale=None):
	'''Standard procedure for turning a possibly even grid into an odd one with a distinct center.'''
	if lenscale == None: lenscale = 1.0
	m,n = shape(dat[0])
	if direct == 0:
		return 1./lenscale*array(dat)[:,slice((1 if m%2==0 else None),None),
			slice((1 if n%2==0 else None),None)]
	elif direct == 1:
		return 1./lenscale*array(dat)[:,slice((-1 if m%2==1 else None),(0 if m%2==0 else None),-1),
			slice((-1 if n%2==1 else None),(0 if n%2==0 else None),-1)]

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	#---fix the height shift
	z0 = 0
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2dgrid(params,x,y):		
	z0,c0,x0,y0,sx,sy,th,poslist = params
	#---fix the height shift
	height = 0
	for pos in poslist:
		height += c0*exp(-((x-x0-pos[0])*cos(th)+(y-y0-pos[1])*sin(th))**2/2./sx**2)*\
			exp(-(-(x-x0-pos[0])*sin(th)+(y-y0-pos[1])*cos(th))**2/2./sy**2)
	return height
		
class ModeCouple():
	'''Class object for holding undulation-curvature coupling data.'''
	def __init__(self):
		self.t2d = [array([]) for i in range(4)]
		self.t2de = [array([]) for i in range(4)]
		self.t1d = [array([]) for i in range(4)]
		self.t1denergy = [array([]) for i in range(4)]
		self.tsum2d = array([])
		self.tsum1d = array([])
		self.tsum1de = array([])
		self.qmagst = array([])
		self.scalefac = 0.0
		self.area = 0.0
		self.c0s = array([])
	def calculate_mode_coupling(self,mset,c0s):
		'''Moved the entire mode-coupling calculation to a single function so it can be repeated easily.'''
		#---compute dimensions
		mset.calculate_undulations(removeavg=removeavg,fitlims=fitlims,
			forcekappa=forcekappa,qmagfilter=qmagfilter)
		grid = mset.griddims
		m,n = grid[0]-1,grid[1]-1
		#---if no curvature field, use zeros
		if c0s == []: c0s = zeros((m,n))
		#---units of hqs are unspecified here but then set by autocorr
		print 'calculating hqs'
		
		#hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
		hqs = []
		for i in range(len(mset.surf)):
			print i
			hqs.append(fftwrap(mset.surf[i])/double(m*n))
		#del mset.surf
		#---units of incoming c0s must be in the units specified by mset.lenscale
		print 'calculating cqs'
		cqs = [fftwrap(mset.lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
		print 2
		Lx,Ly = mean(mset.vecs,axis=0)[0:2]
		print 3
		cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
		print 4
		qmags = mset.lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
			for j in range(0,n)] for i in range(0,m)])
		print 5
		hqsa = autocorr(hqs,direct=0,lenscale=mset.lenscale)
		print 6
		hqsb = autocorr(hqs,direct=1,lenscale=mset.lenscale)
		print 7
		center = [cm,cn]
		print 8
		cqsa = autocorr(cqs,direct=0,lenscale=1.0)
		print 9
		cqsb = autocorr(cqs,direct=1,lenscale=1.0)
		print 10
		qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
		print 11
		mt,nt = shape(hqsa[0])
		print 12
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
		print 13
		scalefac = mset.undulate_kappa*area
		print 14
		#---compute terms which contribute to the elastic energy in 2D
		self.t2d[0] = mean((abs(array(hqsa)))**2,axis=0)
		self.t2d[1] = mean(abs(hqsa)*abs(cqsb),axis=0)
		self.t2d[2] = mean(abs(cqsa*hqsb),axis=0)
		#---VITAL METHODOLOGICAL ELEMENT FOLLOWS
		'''
		The longer you run the mesoscale trajectory, the more "washed out" the curvature field becomes.
		This causes the field to look smaller than the hypothetical fixed curvature field in the CGMD.
		We take a single frame, since it will have the same dimensions.
		You need to think more about the dynamic-ness of the dimple.
		'''
		if 0: self.t2d[3] = mean(abs(cqsa*cqsb),axis=0)
		self.t2d[3] = abs(cqsa*cqsb)[0]
		
		#---compute terms which contribute to the elastic energy in 2D, errors
		self.t2de[0] = std((abs(array(hqsa)))**2,axis=0)
		self.t2de[1] = std(abs(hqsa)*abs(cqsb),axis=0)
		self.t2de[2] = std(abs(cqsa*hqsb),axis=0)
		self.t2de[3] = std(abs(cqsa*cqsb),axis=0)
		#---compute terms which contribute to the elastic energy in 1D
		self.t1d[0] = array([[qmagst[i,j],self.t2d[0][i,j]] for j in range(nt) for i in range(mt)])
		self.t1d[1] = array([[qmagst[i,j],self.t2d[1][i,j]] for j in range(nt) for i in range(mt)])
		self.t1d[2] = array([[qmagst[i,j],self.t2d[2][i,j]] for j in range(nt) for i in range(mt)])
		self.t1d[3] = array([[qmagst[i,j],self.t2d[3][i,j]] for j in range(nt) for i in range(mt)])
		#---sum the terms
		self.tsum2d = scalefac*(self.t2d[0]*qmagst**4-self.t2d[1]*qmagst**2-
			self.t2d[2]*qmagst**2+self.t2d[3])
		self.tsum1d = array([[qmagst[i,j],self.tsum2d[i,j]] for j in range(nt) for i in range(mt)])
		self.tsum2de = array([[sqrt((self.t2de[0][i,j]*qmagst[i,j]**4)**2+
			(qmagst[i,j]**2*self.t2de[1][i,j])**2+
			(qmagst[i,j]**2*self.t2de[2][i,j])**2+
			(self.t2de[3][i,j])**2) 
			for j in range(nt)] for i in range(mt)])
		print 16
		self.tsum1de = array([[qmagst[i,j],
			scalefac*sqrt((self.t2de[0][i,j]*qmagst[i,j]**4)**2+
			(qmagst[i,j]**2*self.t2de[1][i,j])**2+
			(qmagst[i,j]**2*self.t2de[2][i,j])**2+
			(self.t2de[3][i,j])**2)] 
			for j in range(nt) for i in range(mt)])
		#---correct the error for the log scale
		self.tsum1de[self.tsum1de>=self.tsum1d] = self.tsum1d[self.tsum1de>=self.tsum1d]*0.999999
		#---compute terms which are proportional to elastic energy (i.e. scaled by wavevectors)
		self.t1denergy[0] = array([[qmagst[i,j],scalefac*self.t2d[0][i,j]*qmagst[i,j]**4] 
			for j in range(nt) for i in range(mt)])
		self.t1denergy[1] = array([[qmagst[i,j],scalefac*self.t2d[1][i,j]*qmagst[i,j]**2] 
			for j in range(nt) for i in range(mt)])
		self.t1denergy[2] = array([[qmagst[i,j],scalefac*self.t2d[2][i,j]*qmagst[i,j]**2] 
			for j in range(nt) for i in range(mt)])
		self.t1denergy[3] = array([[qmagst[i,j],scalefac*self.t2d[3][i,j]] 
			for j in range(nt) for i in range(mt)])
		#---save variables
		self.area = area
		self.scalefac = scalefac
		self.qmagst = qmagst
		#---store c0s in nm units
		self.c0s = [mset.lenscale*array(c0s[i]) for i in range(len(c0s))]
	
def collapse_spectrum(qs,hs):
	'''For rectangular systems, treat lateral dimensions as equivalent and take their average.'''
	cen = array([i/2 for i in shape(qs)])
	discdists = [[sum(array([i-cen[0],j-cen[1]])**2) 
		for j in range(shape(qs)[1])] for i in range(shape(qs)[0])]
	uniqs = unique(flatten(discdists),return_inverse=True)[1]
	collapsed = [mean(ravel(hs)[where(uniqs==i)]) for i in range(uniqs.max()) 
		if mean(ravel(qs)[where(uniqs==i)]) != 0]
	return collapsed

def equipart_resid(freeparams,both=False):
	'''Define the objective function for matching a simulation to equipartition.'''
	if both:
		bothresids = []
		for n in range(2):
			newkappa,newsigma = freeparams
			mset = msets[n]
			qmagst = mscs[n].qmagst
			area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
				for i in mset.surf_index])/mset.lenscale**2)
			scalefac = newkappa*area
			tsum2d = scalefac*(mscs[n].t2d[0]*qmagst**4-\
				mscs[n].t2d[1]*qmagst**2-mscs[n].t2d[2]*qmagst**2+mscs[n].t2d[3])+\
				1*qmagst**2*mscs[n].t2d[0]*newsigma*area*msets[1].lenscale**2
			xdat = collapse_spectrum(mscs[n].qmagst,mscs[n].qmagst)
			ydat = collapse_spectrum(mscs[n].qmagst,tsum2d)
			bothresids.append([abs(ydat[i]-1)*(1./xdat[i])**1.0 for i in range(len(ydat)) if xdat[i] < 1.0])
		bothresids = array(bothresids)
		return sqrt(bothresids[0]*bothresids[1])
	else:
		newkappa,newsigma = freeparams
		mset = msets[m]
		qmagst = mscs[m].qmagst
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
			for i in mset.surf_index])/mset.lenscale**2)
		scalefac = newkappa*area
		tsum2d = scalefac*(mscs[m].t2d[0]*qmagst**4-\
			mscs[m].t2d[1]*qmagst**2-mscs[m].t2d[2]*qmagst**2+mscs[m].t2d[3])+\
			1*qmagst**2*mscs[m].t2d[0]*newsigma*area*msets[1].lenscale**2
		xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
		ydat = collapse_spectrum(mscs[m].qmagst,tsum2d)
		return [abs(ydat[i]-1)*(1./xdat[i])**1.0 for i in range(len(ydat)) if xdat[i] < 1.0]
		
def construct_hypofield(c0ask):
	#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
	hypo[0] = c0ask
	#---reset the hypothetical C0 field according to the new scaling (replaces the calculation above)
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if analysis_descriptors[a]['simtype'] == 'md':
			anum = analysis_names.index(a)
			mset = msets[anum]
			vecs = mean(mset.vecs,axis=0)
			m,n = mset.griddims
			#---getgrid is xyz points in nm
			getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
				for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
			#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
			params = [0,
				hypo[0]*msets[1].lenscale/mset.lenscale,
				vecs[0]*(1+hypo[1])/mset.lenscale,
				vecs[1]*(1+hypo[2])/mset.lenscale,
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				hypo[5]]
			#---handle curvature fields at the mesoscale which cross the PBC boundary
			c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
				for j in range(3*n)] for i in range(3*m)])
			c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
				for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
			collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
	return [c0hypo for i in range(len(mset.surf))]
	
def bigresid(freeparams):
	newkappa,newsigma,c0ask,extent = freeparams
	hypo = [abs(c0ask),0.5,0.5,extent,extent,0]
	r_2 = extent if extent < 20 else 20
	#---reset the hypothetical C0 field according to the new scaling (replaces the calculation above)
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if analysis_descriptors[a]['simtype'] == 'md':
			anum = analysis_names.index(a)
			mset = msets[anum]
			vecs = mean(mset.vecs,axis=0)
			m,n = mset.griddims
			#---getgrid is xyz points in nm
			getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
				for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
			#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
			params = [0,
				hypo[0]*msets[1].lenscale/mset.lenscale,
				vecs[0]*(1+hypo[1])/mset.lenscale,
				vecs[1]*(1+hypo[2])/mset.lenscale,
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				hypo[5]]
			#---handle curvature fields at the mesoscale which cross the PBC boundary
			c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
				for j in range(3*n)] for i in range(3*m)])
			c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
				for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
			collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
	
	#---local
	mscs = []
	
	#---calculate coupled modes
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		#if 'msc' in globals(): del msc
		msc = ModeCouple()
		msc.calculate_mode_coupling(msets[m],collect_c0s[m])
		mscs.append(msc)
	
	mset = msets[0]
	qmagst = mscs[0].qmagst
	area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
		for i in mset.surf_index])/mset.lenscale**2)
	scalefac = newkappa*area
	tsum2d = scalefac*(mscs[0].t2d[0]*qmagst**4-\
		mscs[0].t2d[1]*qmagst**2-mscs[0].t2d[2]*qmagst**2+mscs[0].t2d[3])+\
		1*qmagst**2*mscs[0].t2d[0]*newsigma*area*msets[1].lenscale**2
	xdat = collapse_spectrum(mscs[0].qmagst,mscs[0].qmagst)
	ydat = collapse_spectrum(mscs[0].qmagst,tsum2d)
	return [abs(ydat[i]-1)*(1./xdat[i])**1.0 for i in range(len(ydat)) if xdat[i] < 1.0]
	
def saddlefunc_span(params,x,y):
	'''Quadrapole saddle function that spans the whole patch.'''
	c0poz,c0neg,sx,sy,x0,y0,dxl,dxr,dyl,dyr = params
	return \
		c0poz*np.exp(-(x-x0-dxl)**2/2/sx**2)+\
		c0poz*np.exp(-(x-x0+dxr)**2/2/sx**2)+\
		-c0neg*np.exp(-(y-y0-dyl)**2/2/sy**2)+\
		-c0neg*np.exp(-(y-y0+dyr)**2/2/sy**2)

def saddlefunc(params,x,y):
	'''Quadrapole saddle function.'''
	c0poz,c0neg,sx,sy,x0,y0,dxx,dyy,flip = params
	fac = -1.0 if flip == True else 1.0
	return \
		-c0neg*np.exp(-(x-x0-dxx)**2/2/sx**2)*\
		np.exp(-(y-y0-fac*dyy)**2/2/sx**2)+\
		c0poz*np.exp(-(x-x0+dxx)**2/2/sx**2)*\
		np.exp(-(y-y0-fac*dyy)**2/2/sx**2)+\
		c0poz*np.exp(-(x-x0-dxx)**2/2/sx**2)*\
		np.exp(-(y-y0+fac*dyy)**2/2/sx**2)+\
		-c0neg*np.exp(-(x-x0+dxx)**2/2/sx**2)*\
		np.exp(-(y-y0+fac*dyy)**2/2/sx**2)
		
#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---connect
if 'conn' not in globals(): 
	try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	except: raise Exception('except: cannot connect to database')
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

#---MANUALLY FITTING DIFFERENT PARAMETERS
#-------------------------------------------------------------------------------------------------------------

if 'compare_parameters_grid' in routine:

	#---settings
	master_details_listing = []
	do_couple_plot = False
	whichtest = ['fullscan','oneplot'][1]
	if whichtest == 'oneplot':
		cgmd_avail_list = ['v616-210000-310000-200','v550-300000-400000-200']
		meso_avail_list = ['v2016-rundir-0','v2016-rundir-4','v2016-rundir-8','v2016-rundir-12']
	elif whichtest == 'fullscan':
		cgmd_avail_list = cgmd_avail
		meso_avail_list = [i for i in analysis_descriptors.keys() if i not in cgmd_avail]

	#---master loop
	master_mscs = []
	for batch_cgmd in cgmd_avail_list:
		for batch_meso in meso_avail_list:

			match_scales = analysis_names = plot_reord = [batch_cgmd,batch_meso]
			status('status: comparing '+batch_cgmd+' '+batch_meso)
			bigname = '-'.join(analysis_names)

			#---BEGIN MATCHING SECTION

			#---allocate if empty
			if 'msets' not in globals(): msets = []
			if 'mscs' not in globals(): mscs = []
			if 'collect_c0s' not in globals(): collect_c0s = []

			#---load and interpolate
			lenscale = 1.0
			for a in analysis_names:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				if 'mset' in globals(): del mset
				mset = MembraneSet()
				if simtype == 'meso':
					params = analysis_descriptors[a]
					#---? type on rundir is int but returns as a float
					pklname = 'pkl.structures.meso.'+\
						params['callsign']+'-'+\
						'C_0-'+str(params['C_0'])+'-'+\
						'len-'+str(params['lenscale'])+'-'+\
						'rundir-'+str(int(params['rundir']))+\
						'.pkl'
					print 'unpickle '+str(pklname)
					mset = unpickle(pickles+pklname)
					print array(mset.surf[0]).max()
					c0s = mset.getdata('c0map').data
					collect_c0s.append(c0s)
					msets.append(mset)
				elif simtype == 'md':
					status('status: loading from MD')
					if 'mset' in globals(): del mset
					mset = unpickle(pickles+locate)
					msets.append(mset)
					#---here we set the hypothetical curvature equal to the 
					#---...induced curvature at the mesoscale
					c0ask = (analysis_descriptors[analysis_names[1]])['C_0']
					hypo[0] = c0ask
					#---compute hypothetical curvature field for the MD simulation
					if hypo != None:
						vecs = mean(mset.vecs,axis=0)
						m,n = mset.griddims
						#---getgrid is xyz points in nm
						getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
							for i in linspace(0,vecs[0]/mset.lenscale,m)])
						#---convert everything to nm
						#---recall z0,c0,x0,y0,sx,sy,th = params
						#---params sets C0 in native units, x0,y0 in proportional units, and sx,sy in nm
						params = [0,hypo[0],
							vecs[0]*hypo[1]/mset.lenscale,
							vecs[1]*hypo[2]/mset.lenscale,
							hypo[3],hypo[4],
							hypo[5]]
						c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)]
							for i in range(m)])
						collect_c0s.append([c0hypo for i in range(len(mset.surf))])
					else:
						collect_c0s.append([])
				elif simtype == 'meso_precomp':
					if 'mset' in globals(): del mset
					mset = unpickle(pickles+locate)
					#---for precomp simulations these units are relative to a0 and already in a reglar grid
					collect_c0s.append(mset.getdata('c0map').data)
					msets.append(mset)

			#---match mesoscale length scales to an MD simulation
			if match_scales != None:
				ref_ind = analysis_names.index(match_scales[0])
				move_ind = analysis_names.index(match_scales[1])
				#---matching the average box vectors here by average and not maximum
				lenscale = mean(mean(msets[move_ind].vecs,axis=0)[:2])/\
					(mean(mean(msets[ref_ind].vecs,axis=0)[:2])/msets[ref_ind].lenscale)
				for a in analysis_names:
					if (analysis_descriptors[a])['simtype'] == 'meso' or \
						(analysis_descriptors[a])['simtype'] == 'meso_precomp':
						#analysis_descriptors[a]['lenscale'] = lenscale
						for i in analysis_descriptors[a]: 
							if i != 'lenscale': vars()[i] = (analysis_descriptors[a])[i]
						msets[analysis_names.index(a)].lenscale = lenscale
				'''
				#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
				hypo[0] = c0ask
				#---reset the hypothetical C0 field according to the new scaling 
				#---... this replaces the calculation above
				for a in analysis_names:
					for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
					if analysis_descriptors[a]['simtype'] == 'md':
						anum = analysis_names.index(a)
						mset = msets[anum]
						vecs = mean(mset.vecs,axis=0)
						m,n = mset.griddims
						#---getgrid is xyz points in nm
						getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
							for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
						#---needs checked, "key step used to have a 0.5*hypo[0] here 
						#---...possibly due to convention"
						params = [0,
							hypo[0]*msets[1].lenscale/mset.lenscale,
							vecs[0]*(1+hypo[1])/mset.lenscale,
							vecs[1]*(1+hypo[2])/mset.lenscale,
							sqrt(r_2)/msets[1].lenscale/sqrt(2),
							sqrt(r_2)/msets[1].lenscale/sqrt(2),
							hypo[5]]
						#---handle curvature fields at the mesoscale which cross the PBC boundary
						c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
							for j in range(3*n)] for i in range(3*m)])
						c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
							for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
						collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
				'''
				construct_hypofield(c0ask)

			#---calculate coupled modes
			for a in analysis_names:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				m = analysis_names.index(a)
				if 'msc' in globals(): del msc
				msc = ModeCouple()
				msc.calculate_mode_coupling(msets[m],collect_c0s[m])
				mscs.append(msc)

			test_kappas = [18,20,22,24,26]
			test_sigmas = [-0.5,0.,0.5,1.0,2.0]
			colors = ['g','b']
			simlabel = ['cgmd','meso']

			#---separate residuals
			if 1:
				#---optimize
				popts = []
				for m in range(2):
					p_opt = leastsq(equipart_resid,array([20.0,0.0]))
					print 'optimized, '+simlabel[m]+' = '+str(list(p_opt[0]))
					test_kappas.append(p_opt[0][0])
					test_sigmas.append(p_opt[0][1])
					popts.append(p_opt)
				master_details_listing.append({
					'batch_meso':batch_meso,
					'batch_cgmd':batch_cgmd,
					'kappa_cgmd':popts[0][0][0],
					'sigma_cgmd':popts[0][0][1],
					'kappa_meso':popts[1][0][0],
					'sigma_meso':popts[1][0][1],
					'resid_meso':mean([i**2 for i in equipart_resid([popts[1][0][0],popts[1][0][1]])]),
					'resid_cgmd':mean([i**2 for i in equipart_resid([popts[0][0][0],popts[0][0][1]])]),
					})
			#---deprecated previously used with equipart_resid(both=True) to do combination residual
			if 0:
				p_opt,cov,infodict,mesg,ie = leastsq(equipart_resid,array([20.0,0.0]),
					full_output=True)
				master_details_listing.append({
					'batch_meso':batch_meso,
					'batch_cgmd':batch_cgmd,
					'kappa':p_opt[0],
					'sigma':p_opt[1],
					'resid':mean([i**2 for i in equipart_resid([p_opt[0],p_opt[1]])]),
					})
			master_mscs.append(mscs)

			#---grid of comparisons
			if whichtest == 'oneplot': 
				gs = gridspec.GridSpec(len(test_kappas),len(test_sigmas),wspace=0.0,hspace=0.0)
				fig = plt.figure(figsize=(12,12))
				for ki in range(len(test_kappas)):
					for si in range(len(test_sigmas)):
						ax = fig.add_subplot(gs[ki,si])
						resid = [0,0]
						for m in range(2):
							newkappa = test_kappas[ki]
							newsigma = test_sigmas[si]
							resid[m] = mean([i**2 for i in equipart_resid([newkappa,newsigma])])
							mset = msets[m]
							qmagst = mscs[m].qmagst
							area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
								for i in mset.surf_index])/mset.lenscale**2)
							scalefac = newkappa*area
							tsum2d = scalefac*(mscs[m].t2d[0]*qmagst**4-\
								mscs[m].t2d[1]*qmagst**2-mscs[m].t2d[2]*qmagst**2+mscs[m].t2d[3])+\
								1*qmagst**2*mscs[m].t2d[0]*newsigma*area*msets[1].lenscale**2
							xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
							ydat = collapse_spectrum(mscs[m].qmagst,tsum2d)
							ax.plot(xdat,ydat,'-',lw=2,label=shortname,c=colors[m],alpha=1.0)
						ax.text(0.05,0.95,
							'$\kappa=$'+'{:.1f}'.format(newkappa)+'\n'+\
							'$\gamma=$'+'{:.1f}'.format(newsigma)+'\n'+\
							'{:.2f}'.format(resid[0])+','+'{:.2f}'.format(resid[1]),
							transform=ax.transAxes,va='top')
						if ki == len(test_kappas)-2 and si == len(test_kappas)-2:
							ax.text(0.05,0.15,'CGMD,opt',transform=ax.transAxes,va='top')
						if ki == len(test_kappas)-1 and si == len(test_kappas)-1:
							ax.text(0.05,0.15,'MESO,opt',transform=ax.transAxes,va='top')
						ax.set_xscale('log')
						ax.set_yscale('log')
						ax.set_ylim((0.1,10))
						ax.set_xlim((0.01,10))			
						ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
						ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
						ax.grid(True)
						ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower',nbins=3))
						ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower',nbins=3))
						if ki < len(test_kappas)-1: ax.set_xticklabels([])
						if si > 0: ax.set_yticklabels([])
				plt.suptitle(analysis_descriptors[batch_cgmd]['label']+' vs '+'meso $C_0='+\
					'{:.3f}'.format(analysis_descriptors[batch_meso]['c_0']*msets[1].lenscale)+\
					'\:{(nm)}^{-1}$',fontsize=fsaxlabel+2)
				plt.savefig(pickles+'fig-coupling-signatures-'+batch_cgmd+'-'+batch_meso+'.png',dpi=300)
				plt.clf()		
			#---clean up the loop
			del mscs,msets

if 'compare_parameters_grid_sweep_plot' in routine:
	magic_lenscale = 0.406905256904
	fig = plt.figure(figsize=(10,6))
	gs = gridspec.GridSpec(1,2)
	c0vkappa = array([[analysis_descriptors[i['batch_meso']]['C_0'],i['kappa_meso']] 
		for i in master_details_listing])
	c0vsigma = array([[analysis_descriptors[i['batch_meso']]['C_0'],i['sigma_meso']] 
		for i in master_details_listing])
	ax = fig.add_subplot(gs[0,0])
	ax.set_title('MESO')
	ax.scatter(magic_lenscale*c0vkappa.T[0],c0vkappa.T[1])
	ax.set_xlabel(r'$C_0\:\mathrm{{nm}^{-1}}$')
	ax.set_ylabel(r'$\kappa$')
	ax = fig.add_subplot(gs[0,1])
	ax.set_title('MESO')
	ax.scatter(magic_lenscale*c0vsigma.T[0],c0vsigma.T[1])
	ax.set_xlabel(r'$C_0\:\mathrm{{nm}^{-1}}$')
	ax.set_ylabel(r'$\sigma$')
	plt.savefig('/home/rpb/example-curvature-kappa-sigma.png')
	plt.clf()

if 'sweep_parameters_cgmd' in routine:

	#---allocate if empty
	if 'msets' not in globals(): msets = []
	if 'mscs' not in globals(): mscs = []
	if 'collect_c0s' not in globals(): collect_c0s = []
	
	analysis_names = [
		'v616-210000-310000-200',
		'v616-10000-209900-200',
		'v550-300000-400000-200',
		][:1]

	#---load and interpolate
	lenscale = 1.0
	for a in analysis_names[:1]:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if 'mset' in globals(): del mset
		mset = MembraneSet()
		if simtype == 'md':
			status('status: loading from MD')
			if 'mset' in globals(): del mset
			mset = unpickle(pickles+locate)
			msets.append(mset)
			#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
			#c0ask = (analysis_descriptors[analysis_names[1]])['C_0']
			c0ask = 0.03
			#'''
			hypo[0] = c0ask
			#---compute hypothetical curvature field for the MD simulation
			if hypo != None:
				vecs = mean(mset.vecs,axis=0)
				m,n = mset.griddims
				#---getgrid is xyz points in nm
				getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
					for i in linspace(0,vecs[0]/mset.lenscale,m)])
				#---convert everything to nm
				#---recall z0,c0,x0,y0,sx,sy,th = params
				#---params sets C0 in native units, x0,y0 in proportional units, and sx,sy in nm
				params = [0,hypo[0],
					vecs[0]*hypo[1]/mset.lenscale,
					vecs[1]*hypo[2]/mset.lenscale,
					hypo[3],hypo[4],
					hypo[5]]
				c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)]
					for i in range(m)])
				collect_c0s.append([c0hypo for i in range(len(mset.surf))])
			else:
				collect_c0s.append([])
			#'''
			#construct_hypofield(c0ask)
			
				
	sweepnames = ['$C_0$','$\kappa$','$\gamma$','$\sigma_{a,b}$']
	solnsets = [
		[[24,-0.3,c0,5] for c0 in arange(0,0.05+0.002,0.002)],
		[[kappa,-0.3,0.02,5] for kappa in arange(15,30+0.5,0.5)],
		[[24,gamma,0.02,5] for gamma in arange(-1.0,4.+0.25,0.25)],
		[[24,-0.3,0.02,extent] for extent in arange(5,40+5,5)],
		]
		
	#---override !!!
	
	solnsets = [
		[[24,-0.3,c0,5] for c0 in [0.02]],
		[[kappa,-0.3,0.02,5] for kappa in arange(15,30+0.5,0.5)],
		[[24,gamma,0.02,5] for gamma in arange(-1.0,4.+0.25,0.25)],
		[[24,-0.3,0.02,extent] for extent in arange(5,40+5,5)],
		]
		
	#solnsets = [
	#	[[24,-0.3,c0,5] for c0 in arange(0,0.05+0.002,0.002)],
	#	]
	
	#---plot
	fig = plt.figure(figsize=(6,10))
	gs = gridspec.GridSpec(4,2)

	#---loop over sweeps
	for sweep in range(len(solnsets)):
		solnset = solnsets[sweep]
		ax = fig.add_subplot(gs[sweep,0])
		ax.set_title(sweepnames[sweep],fontsize=fsaxtitle)
		#---loop over tests
		resids = []
		master_mscs = []
		xdats = []
		ydats = []
		for si in range(len(solnset)):
			soln = solnset[si]
			mscs = []
			newkappa,newsigma,c0ask,extent = soln
			#---reset the hypothetical C0 field according to the new scaling (replaces the calculation above)
			#---note that I replaced msets[1].lenscale with magiclenscale below			
			magic_lenscale = 0.406905256904
			for a in analysis_names[:1]:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				if analysis_descriptors[a]['simtype'] == 'md':
					anum = analysis_names.index(a)
					mset = msets[anum]
					vecs = mean(mset.vecs,axis=0)
					m,n = mset.griddims
					#---getgrid is xyz points in nm
					getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
						for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
					#---set the hypothesis
					hypo = [c0ask,0.5,0.5,extent,extent,0]	
					r_2 = extent
					#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
					params = [0,
						hypo[0]*magic_lenscale/mset.lenscale,
						vecs[0]*(1+hypo[1])/mset.lenscale,
						vecs[1]*(1+hypo[2])/mset.lenscale,
						sqrt(r_2)/magic_lenscale/sqrt(2),
						sqrt(r_2)/magic_lenscale/sqrt(2),
						hypo[5]]
					#---handle curvature fields at the mesoscale which cross the PBC boundary
					c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
						for j in range(3*n)] for i in range(3*m)])
					c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
						for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
					collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
			#---calculate coupled modes
			for a in analysis_names[:1]:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				m = analysis_names.index(a)
				if 'msc' in globals(): del msc
				msc = ModeCouple()
				print 'starting calculate_mode_coupling'
				msc.calculate_mode_coupling(msets[m],collect_c0s[m])
				mscs.append(msc)
			#---compute residuals
			#resid = mean([i**2 for i in equipart_resid([newkappa,newsigma])])
			#resids.append(resid)
			master_mscs.append(mscs)
			mset = msets[m]
			qmagst = mscs[m].qmagst
			area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
				for i in mset.surf_index])/mset.lenscale**2)
			scalefac = newkappa*area
			tsum2d = scalefac*(mscs[m].t2d[0]*qmagst**4-\
				mscs[m].t2d[1]*qmagst**2-mscs[m].t2d[2]*qmagst**2+mscs[m].t2d[3])+\
				1*qmagst**2*mscs[m].t2d[0]*newsigma*area*magic_lenscale**2
			xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat = collapse_spectrum(mscs[m].qmagst,tsum2d)
			xdats.append(xdat)
			ydats.append(ydat)
			ax.plot(xdat,ydat,'-',lw=2,c=mpl.cm.coolwarm(float(si)/len(solnset)))
			pickle.dump([xdat,ydat],open(pickles+'tmp.v616.late.dull.spectrum.pkl','w'))
			#---debug, construction
			raw_input('...')

		#---plot settings
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim((0.1,10))
		ax.set_xlim((0.01,10))			
		ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
		ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
		ax.grid(True)

		#---plot the spectrum differences
		ax = fig.add_subplot(gs[sweep,1])
		for i in range(len(ydats)-1):
			ydiffs = abs(array(ydats[i+1])-array(ydats[i]))
			if sum(ydiffs) > 0:
				thiscolor = mpl.cm.coolwarm(float(i)/(len(ydats)-1)) if sweep == 3 else 'k'
				ax.plot(xdats[i],ydiffs,'-',lw=2,c=thiscolor)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.grid(True)

	#---show
	plt.suptitle(analysis_descriptors[analysis_names[0]]['label']+' parameter sensitivity',fontsize=fsaxlabel+2)
	plt.savefig(pickles+'fig-curvature-coupling-cgmd-sweeps.png')
	plt.show()

#---extra plotting code for memory-adled execution
#ax=plt.subplot(111);ax.plot(xdat,ydat,'-',lw=2,c=mpl.cm.coolwarm(float(si)/len(solnset)));ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim((0.1,10));ax.set_xlim((0.01,10));ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k');ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k');ax.grid(True);plt.show()

#-------------------------------------------------------------------------------------------------------------

def generate_field_saddle(mset,params,ask=None):
	'''Given an mset object and parameters for the saddle function, return the corresponding field.'''
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
		for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
	field_nopbc = array([[saddlefunc(params,getgrid[i,j,0],getgrid[i,j,1]) 
		for j in range(3*n)] for i in range(3*m)])
	field = array(field_nopbc[m:2*m,n:2*n])
	if ask == None: return field
	else: return field/field.max()*ask
	
def generate_field_dimple(mset,params):
	'''Given an mset object and parameters for the saddle function, return the corresponding field.'''
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
		for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
	#---handle curvature fields at the mesoscale which cross the PBC boundary
	field_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
		for j in range(3*n)] for i in range(3*m)])
	field = array([[max([field_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
		for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)])
	return field
	
def generate_field_dimplegrid(mset,params):
	'''Given an mset object and parameters for the saddle function, return the corresponding field.'''
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
		for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
	#---handle curvature fields at the mesoscale which cross the PBC boundary
	field_nopbc = array([[gauss2dgrid(params,getgrid[i,j,0],getgrid[i,j,1]) 
		for j in range(3*n)] for i in range(3*m)])
	field = array([[max([field_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
		for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)])
	return field

def present_field_and_spectrum(mset,params,axl,axr,axr_extras=None,form=None,vmax=None):
	if form == None: raise Exception('except: must specify function type')
	elif form == 'saddle':
		#---set some variables so that you can use verbatim code
		axes = [axl,axr,axr_extras]
		msets = [mset]
		#---prepare hypothesis
		c0ask = params['c_0']/mset.lenscale
		x0,y0 = params['x0'],params['y0']
		wid = params['width']
		hoffset = params['hoffset']
		voffset = params['voffset']
		pozneg = params['pozneg']
		flip = params['flip']
		#---hypothesize
		vecs = mean(msets[0].vecs,axis=0)
		spec_params = [pozneg,1,wid,wid,vecs[0]*(1+x0)/mset.lenscale,vecs[1]*(1+y0)/mset.lenscale,
			hoffset,voffset,flip]
		print flip
		field = generate_field_saddle(mset,spec_params,ask=c0ask)
		fieldsad = field
		#---prepare hypothesis, extra parameters
		kappa,gamma = [24,-0.3]
		#---calculate mode coupling
		msc = ModeCouple()
		msc.calculate_mode_coupling(mset,[field for fr in range(len(mset.surf))])
		#---plot the spectrum
		qmagst = msc.qmagst
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
			for i in mset.surf_index])/mset.lenscale**2)
		scalefac = kappa*area
		tsum2d = scalefac*(msc.t2d[0]*qmagst**4-\
			msc.t2d[1]*qmagst**2-msc.t2d[2]*qmagst**2+msc.t2d[3])+\
			1*qmagst**2*msc.t2d[0]*gamma*area*magic_lenscale**2
		xdat = collapse_spectrum(msc.qmagst,msc.qmagst)
		ydat = collapse_spectrum(msc.qmagst,tsum2d)
		axes[0].plot(xdat,ydat,'-',lw=1.5,c='m',label='saddle')
		if c0ask == 0: c0ask = 0.02
		if vmax == None: vmax = c0ask
		vmax = vmax/mset.lenscale
		axes[1].imshow(field.T,interpolation='nearest',origin='lower',
			vmax=vmax,vmin=-1*vmax,cmap=mpl.cm.RdBu_r)
		if 0:
			axes[1].yaxis.set_label_position("right")
			axes[1].set_ylabel(r'$C_{0,max}='+'{:.3f}'.format(fieldsad.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4,rotation=270)
		else: 
			axes[1].set_title(r'$C_{0,max}='+'{:.3f}'.format(fieldsad.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4)
		if axr_extras != None:
			dat = array(msc.tsum2d)
			cm,cn = [int(i/2) for i in shape(dat)]
			i2wid = 5
			dat = dat[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1]
			im = plotter2d(axes[2],mset,dat=dat,
				lognorm=False,cmap=mpl.cm.jet,inset=False,cmap_washout=1.0,i2wid=i2wid,
				ticklabel_show=[0,0],tickshow=[0,0],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[0.2,2],
				tickskip=int(round(mset.griddims[0]/6,-1)))
	elif form == 'dimplegrid':
		#---set some variables so that you can use verbatim code
		axes = [axl,axr,axr_extras]
		msets = [mset]
		#---prepare hypothesis
		c0ask = params['c_0']/mset.lenscale
		x0,y0 = params['x0'],params['y0']
		extent = params['extent']
		poslist = params['pos']
		hypo = [c0ask,x0,y0,extent,extent,0]
		#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
		vecs = mean(msets[0].vecs,axis=0)
		spec_params = [0,hypo[0],
			vecs[0]*(1+hypo[1])/mset.lenscale,
			vecs[1]*(1+hypo[2])/mset.lenscale,
			sqrt(hypo[3])/magic_lenscale/sqrt(2),
			sqrt(hypo[4])/magic_lenscale/sqrt(2),
			hypo[5],poslist]
		field = generate_field_dimplegrid(mset,spec_params)
		fielddim = field
		#---prepare hypothesis, extra parameters
		kappa,gamma = [24,-0.3]
		#---calculate mode coupling
		msc = ModeCouple()
		msc.calculate_mode_coupling(mset,[field for fr in range(len(mset.surf))])
		#---plot the spectrum
		qmagst = msc.qmagst
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
			for i in mset.surf_index])/mset.lenscale**2)
		scalefac = kappa*area
		tsum2d = scalefac*(msc.t2d[0]*qmagst**4-\
			msc.t2d[1]*qmagst**2-msc.t2d[2]*qmagst**2+msc.t2d[3])+\
			1*qmagst**2*msc.t2d[0]*gamma*area*magic_lenscale**2
		xdat = collapse_spectrum(msc.qmagst,msc.qmagst)
		ydat = collapse_spectrum(msc.qmagst,tsum2d)
		axes[0].plot(xdat,ydat,'-',lw=1.5,c='m',label='dimple')
		if c0ask == 0: c0ask = 0.02
		if vmax == None: vmax = c0ask
		vmax = vmax/mset.lenscale
		axes[1].imshow(field.T,interpolation='nearest',origin='lower',
			vmax=vmax,vmin=-1*vmax,cmap=mpl.cm.RdBu_r)
		if 0:
			axes[1].yaxis.set_label_position("right")
			axes[1].set_ylabel(r'$C_{0,max}='+'{:.3f}'.format(fielddim.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4,rotation=270)
		else: 
			axes[1].set_title(r'$C_{0,max}='+'{:.3f}'.format(fielddim.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4)
		if axr_extras != None:
			dat = array(msc.tsum2d)
			cm,cn = [int(i/2) for i in shape(dat)]
			i2wid = 5
			dat = dat[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1]
			im = plotter2d(axes[2],mset,dat=dat,
				lognorm=False,cmap=mpl.cm.jet,inset=False,cmap_washout=1.0,i2wid=i2wid,
				ticklabel_show=[0,0],tickshow=[0,0],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[0.2,2],
				tickskip=int(round(mset.griddims[0]/6,-1)))
	elif form == 'dimple':
		#---set some variables so that you can use verbatim code
		axes = [axl,axr,axr_extras]
		msets = [mset]
		#---prepare hypothesis
		c0ask = params['c_0']/mset.lenscale
		x0,y0 = params['x0'],params['y0']
		extent = params['extent']
		hypo = [c0ask,x0,y0,extent,extent,0]
		#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
		vecs = mean(msets[0].vecs,axis=0)
		spec_params = [0,hypo[0],
			vecs[0]*(1+hypo[1])/mset.lenscale,
			vecs[1]*(1+hypo[2])/mset.lenscale,
			sqrt(hypo[3])/magic_lenscale/sqrt(2),
			sqrt(hypo[4])/magic_lenscale/sqrt(2),
			hypo[5]]	
		field = generate_field_dimple(mset,spec_params)
		fielddim = field
		#---prepare hypothesis, extra parameters
		kappa,gamma = [24,-0.3]
		#---calculate mode coupling
		msc = ModeCouple()
		msc.calculate_mode_coupling(mset,[field for fr in range(len(mset.surf))])
		#---plot the spectrum
		qmagst = msc.qmagst
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
			for i in mset.surf_index])/mset.lenscale**2)
		scalefac = kappa*area
		tsum2d = scalefac*(msc.t2d[0]*qmagst**4-\
			msc.t2d[1]*qmagst**2-msc.t2d[2]*qmagst**2+msc.t2d[3])+\
			1*qmagst**2*msc.t2d[0]*gamma*area*magic_lenscale**2
		xdat = collapse_spectrum(msc.qmagst,msc.qmagst)
		ydat = collapse_spectrum(msc.qmagst,tsum2d)
		axes[0].plot(xdat,ydat,'-',lw=1.5,c='g',label='dimple')
		if c0ask == 0: c0ask = 0.02
		if vmax == None: vmax = c0ask
		vmax = vmax/mset.lenscale
		axes[1].imshow(field.T,interpolation='nearest',origin='lower',
			vmax=vmax,vmin=-1*vmax,cmap=mpl.cm.RdBu_r)
		if 0:
			axes[1].yaxis.set_label_position("right")
			axes[1].set_ylabel(r'$C_{0,max}='+'{:.3f}'.format(fielddim.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4,rotation=270)
		else: 
			axes[1].set_title(r'$C_{0,max}='+'{:.3f}'.format(fielddim.max()*mset.lenscale)+'$',
				fontsize=fsaxlabel-4)
		if axr_extras != None:
			dat = array(msc.tsum2d)
			cm,cn = [int(i/2) for i in shape(dat)]
			i2wid = 5
			dat = dat[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1]
			im = plotter2d(axes[2],mset,dat=dat,
				lognorm=False,cmap=mpl.cm.jet,inset=False,cmap_washout=1.0,i2wid=i2wid,
				ticklabel_show=[0,0],tickshow=[0,0],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[0.2,2],
				tickskip=int(round(mset.griddims[0]/6,-1)))

	
#---analyze the CGMD simulations with different hypothetical curvature fields
if 'cgmd_test_newfield' in routine:

	#---temporary setting
	magic_lenscale = 0.406905256904
	
	#---additional 2D spectrum plots
	show2dspec = True

	#---settings
	analysis_names = [
		'v616-210000-310000-200',
		'v550-300000-400000-200',
		][:1]	

	#---load
	msets = []
	lenscale = 1.0
	for a in analysis_names[:1]:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if 'mset' in globals(): del mset
		mset = MembraneSet()
		mset = unpickle(pickles+locate)
		msets.append(mset)
	mset = msets[0]
		
	#---prepare sweep
	basic_saddle = {
		'c_0':0.0,
		'hoffset':5,
		'voffset':5,
		'x0':0.5,'y0':0.5,
		'width':5,
		'pozneg':2,
		'flip':False}
	sweep_saddle = []
	curvlist = arange(0,0.05,0.01)
	for i in range(len(curvlist)):
		sweep_saddle.append(dict(basic_saddle))
		sweep_saddle[-1]['c_0'] = curvlist[i]
		sweep_saddle[-1]['flip'] = True if i%2 == 0 else False

	#---prepare sweep
	basic_dimple = {
		'c_0':0.0,
		'extent':30,
		'x0':0.5,'y0':0.5,}
	sweep_dimple = []
	curvlist = arange(0,0.05,0.01)
	for i in range(len(curvlist)):
		sweep_dimple.append(dict(basic_dimple))
		sweep_dimple[-1]['c_0'] = curvlist[i]
		
	#---prepare sweep
	dd = 10.
	th = -sqrt(3)/2*dd
	basic_dimplegrid = {
		'c_0':0.0,
		'extent':3,
		'x0':0.5,'y0':0.5,
		'pos':
		[(-1*dd,0),(0,0),(dd,0),(-dd/2,th),(dd/2,th),(dd*3/2,th),(-dd/2,-th),(dd/2,-th)]
		}
	sweep_dimplegrid = []
	curvlist = arange(0,0.05,0.01)
	for i in range(len(curvlist)):
		sweep_dimplegrid.append(dict(basic_dimplegrid))
		sweep_dimplegrid[-1]['c_0'] = curvlist[i]
	
	#---prepare figure
	fig = plt.figure(figsize=(10,3*len(sweep_saddle)))
	ncols = 5 if show2dspec else 3
	gs = gridspec.GridSpec(len(sweep_saddle),ncols)
	axes = [[fig.add_subplot(gs[j,i]) for i in range(ncols)] for j in range(len(sweep_saddle))]
	
	for i in range(len(sweep_saddle)):
		if 0:
			present_field_and_spectrum(mset,sweep_saddle[i],
				axes[i][0],axes[i][1],form='saddle',vmax=0.05,
				axr_extras=(axes[i][3] if show2dspec else None))
		else:
			present_field_and_spectrum(mset,sweep_dimplegrid[i],
				axes[i][0],axes[i][1],form='dimplegrid',vmax=0.05,
				axr_extras=(axes[i][3] if show2dspec else None))
		present_field_and_spectrum(mset,sweep_dimple[i],
			axes[i][0],axes[i][2],form='dimple',vmax=0.05,
			axr_extras=(axes[i][4] if show2dspec else None))

		#---plot settings
		ax = axes[i][0]
		ax.set_xlim((0.01,10))
		ax.set_ylim((0.1,10))
		ax.set_aspect(3./2.)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
		ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
		ax.grid(True)
		
	plt.suptitle('saddle vs dimple, '+\
		analysis_descriptors[analysis_names[0]]['label'],fontsize=fsaxlabel+2)
	plt.savefig(pickles+'fig-saddle_v_dimple-'+analysis_names[0]+'.png',dpi=300)
	plt.show()


#---jot-sharp.py code follows
if 0:
	import time

	if 0:
		self = mset
		frameno = 10
		rounder = 1
		interp='best'
		self.gotoframe(frameno)
		topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
			for i in self.monolayer_residues[0]])
		topxyz = topxyz - mean(topxyz)
		if 1: topxyzwrap = self.wrappbc(topxyz,vecs=self.vec(frameno),mode='grow')
		else: topxyzwrap = topxyz
		self.griddims = [int(round(self.vec(1)[0]/rounder)),int(round(self.vec(1)[1]/rounder))]
		st = time.time()
		topmesh = self.makemesh(topxyzwrap,self.vec(frameno),self.griddims,method=interp)
		print 1./60*(time.time()-st)
		tmp = mset.rezipgrid(topmesh)
		print 1./60*(time.time()-st)
		tmp = tmp-mean(tmp)
		if 0: plt.imshow(array(tmp).T,origin='lower',interpolation='nearest');plt.show()
		st = time.time()
		fineness = 500j,500j
		gridx,gridy = mgrid[0:self.vec(frameno)[0]:fineness[0],0:self.vec(frameno)[1]:fineness[1]]
		points = topmesh[:,:2]
		values = topmesh[:,2]
		gridz2 = scipy.interpolate.griddata(points, values, (gridx, gridy), method='cubic')
		print 1./60*(time.time()-st)
		plt.imshow(tile(array(gridz2).T,(3,3)),origin='lower',interpolation='nearest',cmap=mpl.cm.binary);plt.show()

	if 0:
		newway = True

		self = mset
		frameno = 10
		rounder = 1
		interp='sharp_cubic_spline'
		self.gotoframe(frameno)
		topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
			for i in self.monolayer_residues[0]])
		topxyz = topxyz - mean(topxyz)
		if 1: topxyzwrap = self.wrappbc(topxyz,vecs=self.vec(frameno),mode='grow')
		else: topxyzwrap = topxyz
		self.griddims = [int(round(self.vec(1)[0]/rounder)),int(round(self.vec(1)[1]/rounder))]
		st = time.time()
		if newway:
			topmesh = self.makemesh(topxyzwrap,self.vec(frameno),self.griddims,method=interp,fine=1000)
		else:
			topmesh = self.makemesh(topxyzwrap,self.vec(frameno),self.griddims,method='best',fine=1000)
		print 1./60*(time.time()-st)
		tmp = mset.rezipgrid(topmesh)
		tmp = tmp-mean(tmp)
		if newway:	
			st = time.time()
			fineness = 500j,500j
			gridx,gridy = mgrid[0:self.vec(frameno)[0]:fineness[0],0:self.vec(frameno)[1]:fineness[1]]
			points = topxyzwrap[:,:2]
			values = topxyzwrap[:,2]
			gridz2 = scipy.interpolate.griddata(points, values, (gridx, gridy), method='cubic')
			print 1./60*(time.time()-st)
			plt.imshow(tile(array(gridz2).T,(3,3)),origin='lower',interpolation='nearest',cmap=mpl.cm.binary);plt.show()
		else: 
			plt.imshow(tile(array(tmp).T,(3,3)),origin='lower',interpolation='nearest',cmap=mpl.cm.binary);plt.show()
		

	if 0:
		from membrainrunner import *
		execfile('locations.py')
		mset = unpickle(pickles+'pkl.structures.sharp.space2A.membrane-v616.u6-trestles.md.part0004.110000-209900-200.pkl')
		mset.calculate_undulations()
		plotter_undulate(mset);plt.show()
		plt.imshow(mset.surf[0]);plt.show()



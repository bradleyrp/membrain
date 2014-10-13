#!/usr/bin/python -i

from membrainrunner import *
allsets = dict()
execfile('locations.py',allsets)
execfile('locations/header.py',allsets)
from module_modecouple.ModeCouple import *
from module_modecouple.ModeSum import *
from module_database.DataFace import DataFace
import copy
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---choose what to do
routine = [
	'calc',
	'hypothesize',
	'coupling_solution',
	][1:]

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
		'C_0':0.000,
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
		'C_0':[0.000,0.001,0.005,0.01,0.018,0.02,0.022,0.024,0.03,0.035,0.04,0.05]
		},
	}

#---sweep over bending rigidity
sweep_kappa = {
	'kappa':{
		'fore':[20,22,24,28,32],
		},
	}

#---HYPOTHESIS
#-------------------------------------------------------------------------------------------------------------

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

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def specfilter_calc(specs):
	''' '''
	filt = dict()
	for key in specs.keys():
		if key in hypothesis_calc.keys():
			filt[key] = specs[key]
	return filt

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def plot_spectrum(dat,denote,logplot=False,show=False):
	'''Temporary spectrum plotting program.'''
	datflat = abs(array(dat))+(0.*10**-10 if logplot else 0)
	print datflat.max()
	print datflat.min()
	ax = plt.subplot(111)
	im = ax.imshow(
		datflat.T,
		interpolation='nearest',
		origin='lower',
		norm=(mpl.colors.LogNorm() if logplot else None))
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	if not show: 
		plt.savefig(
			allsets['pickles']+'FIGURES/fig-modecouple-'+denote+\
			('-log-1' if logplot else '')+'.png',dpi=300)
		plt.clf()
	else: plt.show()
	
def plot_spectrum1d(xvals,yvals,denote,show=False,logplot=True,ylims=None):
	'''Temporary spectrum plotting program.'''
	ax = plt.subplot(111)
	ax.plot(xvals,yvals,'o')
	ax.set_xscale('log')
	ax.set_yscale('log')	
	if ylims != None: 
		ax.set_ylim(ylims)
	#---fit line
	color = 'r'
	label = ''
	kappa = 20
	qmagfilter = [0.15,0.5]
	specfilter = array(filter(lambda x: x[0] >= qmagfilter[0] and x[0] 
		<= qmagfilter[1],array([xvals,yvals]).T))
	leftcom = [mean(log(specfilter[:,0])),mean(log(specfilter[:,1]))]
	[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
	az_enforced = leftcom[1]+2.*leftcom[0]
	ax.plot(specfilter[:,0],[exp(az)*(i**bz) 
		for i in specfilter[:,0]],c=color,lw=1.5,
		label=(None if label == None else label+'\n'+\
		r'$\boldsymbol{\kappa} = '+str('%3.1f'%kappa)+'\:k_BT$'))
	print [bz,az]
	#---save
	if not show: 
		plt.savefig(
			allsets['pickles']+'FIGURES/fig-modecouple-'+denote+\
			('-log-1' if logplot else '')+'.png',dpi=300)
		plt.clf()
	else: plt.show()
	
def perfect_collapser(xs,ys):
	'''Return two arrays for a "perfect" collapse.'''
	xsort = array(sort(list(set(xs)))[1:])
	inds = argmax(array([xs==xsort[i] for i in range(len(xsort))]),axis=0)
	col = array([mean(ys[where(inds==i)]) for i in range(len(xsort))])
	return xsort,col
	
def blurry_binner(xs,ys,binwidth=0.05):
	'''Group wavevectors by bins.'''
	blurred = (xs/binwidth).astype(int)
	xsort = array(sort(list(set(blurred)))[1:])
	inds = argmax(array([(xs/binwidth).astype(int)==xsort[i] for i in range(len(xsort))]),axis=0)
	coly = array([mean(ys[where(inds==i)]) for i in range(len(xsort))])
	colx = array([mean(xs[where(inds==i)]) for i in range(len(xsort))])
	return colx,coly,inds
	

#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if 'calc' in routine:

	'''
	Since parts of the hypothesis, such as the bending rigidity field, are easy to compute, the calculation
	section doesn't require a full hypothesis. Instead, it computes the terms which contribute to the final
	Helfrich sum and depend on a minimal hypothesis, in this case consisting only of the curvature field. 
	These results are then stored in hd5f binaries which are unpacked during the "hypothesize" stage.
	'''
	
	for hypothesis in hypotheses:
	
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

if 'hypothesize' in routine and 'ms' not in globals():

	#---select a hypothesis via the pass-through argument to the script
	hypothesis = hypotheses[int(args.passed)]

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
		c0qs_c0qps=termlist[3],
		do_diag=False)
	
	#---check on the kappa field	
	if 0: plt.imshow(real(ms.kqs).T,interpolation='nearest',origin='lower');plt.show()
	
elif 'ms' in globals(): print 'refusing to continue because ms already defined'
		
#---TESTING
#-------------------------------------------------------------------------------------------------------------

if 'testing_basic_undulation_spectra' in routine:

	#---redoing termlist from scratch
	if 'termdat_test' not in globals():

		m2,n2 = shape(ms.hqs)[1:]-0*array([1,1])
		st = time.time()
		termdat_test = zeros((m2*n2,m2*n2))
		for fr in range(len(ms.hqs)):
			status('fr = '+str(fr),start=st,i=fr,looplen=len(ms.hqs))
			#---note also removed trim from here
			#---note that this tests the outer function
			termdat_test += abs(outer(ms.hqs[fr],ms.hqs[fr]))/float(len(ms.hqs))
			if 0: termdat_test += outer(abs(ms.hqs[fr]),abs(ms.hqs[fr]))/float(len(ms.hqs))
		
	#---remake of the built-in undulation calculator
	#---...which perfectly captures the original code and correct scaling
	#---...without cutting any of the original heights from the edges of the grid
	trim = 0
	m,n = m2-trim,n2-trim
	lxy = ms.mset.vecs
	yv,xv = meshgrid(range(n),range(m))
	undulate_raw = [fft.fftshift(fft.fft2(
		ms.mset.surf[fr][:-1 if trim == 1 else None,:-1 if trim == 1 else None])) 
		for fr in range(len(ms.mset.surf))]
	undulate_hqhq2d = mean(array(1.*(abs(array(undulate_raw))/double(m*n))**2),axis=0)
	spec_center = [int(i/2) for i in shape(undulate_raw)[1:]]
	qs = ms.mset.lenscale*array([(sqrt((2*pi*(array(xv)-spec_center[0])/lxy[f][0])**2+
		(2*pi*(array(yv)-spec_center[1])/lxy[f][1])**2)) for f in range(len(lxy))])
	undulate_qmag2d = mean(qs,axis=0)
	ys = reshape(undulate_hqhq2d,-1)[1:]
	xs = reshape(undulate_qmag2d,-1)[1:]
	print 'plotting the undulations via the original method...'
	plot_spectrum1d(xs,ys,'',show=True,logplot=True,ylims=[10**-6,10**2])	
	
	#---new version from outer to confirm that outer works to compute the fluctuations
	Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]		
	qmags = ms.mset.lenscale*array([[sqrt(((i)/((Lx)/1.)*2*pi)**2+((j)/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])	
	qvals = reshape(qmags,-1)
	qmagsshift = ms.mset.lenscale*array([[sqrt(((i-m2*(i>m2/2))/((Lx)/1.)*2*pi)**2+
		((j-n2*(j>n2/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])
	xvals = reshape(qmagsshift,-1)[1:]
	#---accidentally deleted the definition of tmp a couple weeks ago
	#---replaced it with an obvious choice here but this probably needs checked
	tmpa = outer(ms.hqs[0],ms.hqs[0])
	tmp = [abs(tmpa[i,i]) for i in range(len(qvals))]
	#---corrected method that uses the whole termdat_test
	tmp = array([termdat_test[i,i] for i in range(len(qvals))])
	yvals = tmp[1:]
	#---plot the spectrum after collapsing in three different ways	
	#---raw data
	print 'plotting raw data...'
	plot_spectrum1d(xvals,yvals,'',show=True,logplot=True,)	
	#---perfect binner
	print 'plotting data after perfect binning...'
	xvals2,yvals2 = perfect_collapser(xvals,yvals)
	plot_spectrum1d(xvals2,yvals2,'',show=True,logplot=True,)	
	#---blurry binner
	print 'plotting data after blurry binning...'
	xvals3,yvals3,inds = blurry_binner(xvals,yvals)
	plot_spectrum1d(xvals3,yvals3,'',show=True,logplot=True,)	
		
#---SOLUTION
#-------------------------------------------------------------------------------------------------------------

if 'coupling_solution' in routine:

	#---collect the full matrix from the ModeSum object
	termdat_test = ms.full_matrix

	#---prepare the diagonalto set the bin sizes
	m2,n2 = shape(ms.hqs)[1:]-0*array([1,1])
	Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]		
	qmags = ms.mset.lenscale*array([[sqrt(((i)/((Lx)/1.)*2*pi)**2+((j)/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])	
	qvals = reshape(qmags,-1)
	qmagsshift = ms.mset.lenscale*array([[sqrt(((i-m2*(i>m2/2))/((Lx)/1.)*2*pi)**2+
		((j-n2*(j>n2/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])
	xvals = reshape(qmagsshift,-1)[1:]
	diag = array([termdat_test[i,i] for i in range(len(qvals))])
	yvals = diag[1:]
	xvals3,yvals3,inds = blurry_binner(xvals,yvals)

	#---reduce dimensionality and average
	reddim = inds.max()
	tiny = zeros((reddim,reddim))
	tinycount = zeros((reddim,reddim))
	st = time.time()
	for i in range(reddim):
		status('i = '+str(i),start=st,i=i,looplen=reddim)
		for j in range(reddim):
			#---I did one check, but this needs another to make sure indexing is working right
			jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
			sieve = termdat_test[jjj[...,0],jjj[...,1]]
			tiny[i,j] = mean(sieve)
			tinycount[i,j]  = len(sieve)
	plt.imshow(array(tiny/tinycount).T,origin='lower',interpolation='nearest',norm=mpl.colors.LogNorm())
	plt.savefig(allsets['pickles']+'COUPLESOLN-'+str(args.passed)+'.png')
	plt.show()
	



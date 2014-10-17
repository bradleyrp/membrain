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
	'testing_basic_undulation_spectra',
	'coupling_solution_loops',
	][-1:]

#---select a simulation system
callsign = [
	'v550-300000-400000-200',
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	][2]
	
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

#---sweep is an arbitrary list of couples for doing multidimensional parameter sweeps
#---...the first element of the couple is the "route" to the final value (i.e. sequential dict keys)
#---...and the second element of the couple is the list of parameters to sweep over
#---...all combinations will be combined into a master hypothesis list
sweep_calc = [
	[['curvature','C_0'],[0.000,0.001,0.005,0.01,0.018,0.02,0.022,0.024,0.03,0.035,0.04,0.05]],
	[['kappa','fore'],[20,22,24,28,32]],
	[['kappa','back'],[18,20]],
	][:2]
	
#---sweep for analysis
sweep_look = [
	[['curvature','C_0'],[0.000,0.001,0.005,0.01,0.018,0.02,0.022,0.024,0.03,0.035,0.04,0.05]],
	[['kappa','fore'],[20,22,24,28,32]],
	[['kappa','radius'],[5,10,15,20,25]],	
	]
	
#---set the binwidth for the "blurry" binner function
bin_width = 0.1

#---HYPOTHESIS
#-------------------------------------------------------------------------------------------------------------

sweep = sweep_calc if 'calc' in routine else sweep_look

#---extract a list of lists of parameters to sweep over
t = [i[1] for i in sweep]
#---take all combinations via meshgrid
allcombos = reshape(array(meshgrid(*t)).T,(-1,3))

#---assemble a list of hypotheses from all possible combinations of the sweep values
#---note that this code is general, and works for an arbitrarily deep dictionary
hypotheses = []
#---for each combo generate a new hypothesis
for combo in allcombos:
	#---start with the default hypothesis
	newhypo = copy.deepcopy(hypothesis_default)
	#---each combo has a value and a route which is a sequence of dictionary keys
	#---...we loop over each route to set each final value for the sweep
	for routenum in range(len(sweep)):
		#---to get to the deepest part of that route we use tmp as a pointer
		#---...and iteratively traverse one level until the second to last level
		tmp = newhypo[sweep[routenum][0][0]]
		for i in sweep[routenum][0][1:-1]: tmp = tmp[i]
		#---at the final level, we now have a pointer to the lowest dictionary where we set the value
		tmp[sweep[routenum][0][-1]] = combo[routenum]
	#---once we set all the values, the hypothesis is ready
	hypotheses.append(newhypo)
		
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
	
def blurry_binner(xs,ys,bin_width=0.05):
	'''Group wavevectors by bins.'''
	blurred = (xs/bin_width).astype(int)
	xsort = array(sort(list(set(blurred)))[1:])
	inds = argmax(array([(xs/bin_width).astype(int)==xsort[i] for i in range(len(xsort))]),axis=0)
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
	if 0: plt.imshow(real(ms.kfield).T,interpolation='nearest',origin='lower');plt.show()
	
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
	xvals3,yvals3,inds = blurry_binner(xvals,yvals,bin_width=bin_width)
	plot_spectrum1d(xvals3,yvals3,'',show=True,logplot=True,)	
		
#---SOLUTION
#-------------------------------------------------------------------------------------------------------------

if 'coupling_solution' in routine:

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
	
	#---"debugging" the math here
	qs = array(meshgrid(arange(m2),arange(n2))).T
	lenscale = ms.mset.lenscale
	us = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
	qmagstd = sqrt(sum((us/array([Lx,Ly])*2*pi*lenscale)**2,axis=1))
	#---note that qdq is equal to ms.qdotq
	qdq = outer(qmagstd,qmagstd)
	if not all(qdq == ms.qdotq): raise Exception('except: problem with wavevectors')
	qd = [qdq[i,i] for i in range(len(qdq))]
	#---recap from earlier
	qmagsshift = ms.mset.lenscale*array([[sqrt(((i-m2*(i>m2/2))/((Lx)/1.)*2*pi)**2+
		((j-n2*(j>n2/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])
	xvals = reshape(qmagsshift,-1)[1:]
	#---redefine
	qdq = outer(qmagsshift,qmagsshift)
	qdqd = qdq[range(len(qdq)),range(len(qdq))]


	#---collect the full matrix from the ModeSum object
	if 0: termdat_test = ms.full_matrix
	#---"debugging" the math here by fixing qdotq	
	termdat_test = (termlist[0]*qdq*qdq+termlist[1]*qdq+termlist[2]+termlist[3])*abs(ms.kqqp)

	diag = array([termdat_test[i,i] for i in range(len(qvals))])
	yvals = diag[1:]
	
	#---having used the diagonal of termdat_test to establish the "simple" fluctuation spectrum
	#---...we now use the blurry binner to retrieve the wavevector magnitudes after reducing dimensionality
	xvals3,yvals3,inds = blurry_binner(xvals,yvals,bin_width=bin_width)

	#---the maximum index returned by blurry_binner is equal to the number of reduced dimensions
	reddim = inds.max()
	
	if 'tiny' not in globals():

		#---tiny averages all magnitudes of the target data (termdat_test) with wavevectors in the i,j bin
		tiny = zeros((reddim,reddim))
		#---tinycount records the number of values that contributed to the mean
		tinycount = zeros((reddim,reddim))
	
		#---loop over bins and take average
		st = time.time()
		for i in range(reddim):
			status('i = '+str(i),start=st,i=i,looplen=reddim)
			for j in range(reddim):
				#---check code once
				if i == 5 and j == 3:
					#---the jjj matrix provides a matrix of combinations of positions in inds equal to i,j
					jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
					#---jj provides another way to extract the data from termdat_test via a 1D array
					jj = reshape(jjj,(product(shape(jjj)[:2]),2))
					#---we can use the 1D version and reshape it to get the "sieve" of correct values
					sieve_via_reshape = reshape(termdat_test[jj[:,0],jj[:,1]],shape(jjj)[:2])
					#---but the 1D method is redundant because we can access it directly as well
					sieve_via_direct = termdat_test[jjj[...,0],jjj[...,1]]
					#---we check that both methods agree
					reshape_passes = all(sieve_via_reshape==sieve_via_direct)
					#---and more importantly, we confirm that e.g. the last element in jjj is correct
					example_passes = all([inds[jjj[-1,-1][zz]]==[i,j][zz] for zz in range(2)])
					#---break if the test fails
					if not reshape_passes or not example_passes: raise Exception('except: code broke')				
				#---find indices i,j in the inds list and take all combinations via transposed meshgrid
				jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
				#---here we "sieve" the correct values from termdat_test
				sieve = termdat_test[jjj[...,0],jjj[...,1]]
				#---average the termdat_test values with wavevectors with indices i,j
				tiny[i,j] = mean(sieve)
				tinycount[i,j] = len(sieve)

		#---compute eigenspectrum
		fullans = linalg.eig(abs(tiny))
	
	import matplotlib.gridspec as gridspec

	fig = plt.figure(figsize=(16,8))
	gsij = 4,3
	gs = gridspec.GridSpec(gsij[1],gsij[0])
	axes = [[plt.subplot(gs[i,j]) if 1 else None 
		for j in range(gsij[0])] for i in range(gsij[1])]
	if 0: plt.suptitle('Coupling results, 8xENTH')

	#---full matrix before eigenspectrum
	ax = axes[0][0]
	im = ax.imshow(array(tiny).T,origin='lower',
		interpolation='nearest',norm=mpl.colors.LogNorm())
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	ax.set_title('reduced dimensions')
	ax.set_xlabel('bin number')
	ax.set_ylabel('bin number')
		
	#---eigenvectors
	ax = axes[0][1]
	im = ax.imshow(array(fullans[1]).T,origin='lower',interpolation='nearest',
		vmax=0.8,vmin=-0.8,cmap=mpl.cm.RdBu)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	ax.set_title('eigenvectors')
	ax.set_ylabel('eigenvector number')
	ax.set_xlabel('bin number')

	#---eigenspectrum
	ax = axes[1][0]
	ax.plot(range(len(fullans[0])),abs(fullans[0]))
	ax.set_title('eigenspectrum')
	ax.set_yscale('log')
	ax.set_xlabel('eigenvector number')
	ax.grid(True)

	#---top eigenvalues
	ax = axes[1][1]
	opacs = [0.2,0.3,0.4,0.5,1]
	tops = [1,3,10,20,len(fullans[1])]
	prefactor = 2*product(mean(ms.mset.vecs,axis=0)[:2]/10.)
	for m in range(len(tops)):
		many = tops[m]
		combo = mean([fullans[1][i]*fullans[0][i] for i in range(many)],axis=0)
		ax.plot(xvals3[1:],prefactor*abs(combo),alpha=opacs[m])
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('wavevector ${nm}^{-1}$')
	ax.set_title('summed eigenvectors')
	ax.grid(True)

	#---kappa field
	if 0:
		ax = axes[1][2]
		im = ax.imshow(real(ms.kfield).T,interpolation='nearest',origin='lower')
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		ax.set_title('kappa field')
		ax.set_xlabel('x')
		ax.set_ylabel('y')
	#---kappa field FFT
	else:
		ax = axes[1][2]
		im = ax.imshow(real(ms.kqs).T,interpolation='nearest',origin='lower',
			norm=mpl.colors.LogNorm())
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		ax.set_title('kappa field FFT')
		ax.set_xlabel('x')
		ax.set_ylabel('y')

	#---curvature field
	ax = axes[0][2]
	im = ax.imshow(real(ms.hfield).T,interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	ax.set_title('curvature field')
	ax.set_xlabel('x')
	ax.set_ylabel('y')

	#---term 0
	ax = axes[2][0]
	im = ax.imshow(array(termlist[0]).T,origin='lower',interpolation='nearest',
		norm=mpl.colors.LogNorm())
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	ax.set_title('hqhq raw term')
	
	#---full_matrix
	ax = axes[2][2]
	im = ax.imshow(array(ms.full_matrix).T,origin='lower',interpolation='nearest',
		norm=mpl.colors.LogNorm(),vmin=10**-5,vmax=10**5)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	ax.set_title('full matrix')

	#---term 1
	if hypothesis['curvature']['C_0'] != 0:
		ax = axes[2][1]
		im = ax.imshow(array(termlist[1]).T,origin='lower',interpolation='nearest',
			norm=mpl.colors.LogNorm())
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		ax.set_title('cqhq raw term')
	else: fig.delaxes(axes[2][1])

	#---previous tests to show that you get q4 scaling regardless of the order
	#---...and also to debug a problem with the construction of the q values
	if 0:	
		ax = axes[0][3]
		diag = array([termlist[0][i,i] for i in range(len(qvals))])
		ax.plot(xvals,diag[1:],'b.')
		diag = array([termlist[0][i,i] for i in range(len(qvals))])
		ax.plot(xvals,diag[1:]*xvals**4,'r.')
		ax.set_title('hqhq diagonal')
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_xlabel('wavevector ${nm}^{-1}$')
		ax.grid(True)	
	
		ax = axes[1][3]
		tmp = qdqd*qdqd*termlist[0]
		diag = array([tmp[i,i] for i in range(len(qvals))])
		ax.plot(xvals,diag[1:],'b.')
		ax.set_title(r'$\mathrm{({hqhq}{q}^{4})_{diag}}$')
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_xlabel('wavevector ${nm}^{-1}$')
		ax.grid(True)
	else:
		fig.delaxes(axes[1][3])
		fig.delaxes(axes[2][3])
		
	#---look at how the correction terms affect the final result
	ax = axes[0][3]
	diag = [prefactor*array([termlist[j][i,i] for i in range(len(qvals))]) for j in range(4)]
	kdiag = array([ms.kqqp[i,i] for i in range(len(qvals))])
	clrs = ['r','g','g','b']
	facs = [xvals**4,xvals**2,xvals**2,1.]
	for k in range(4): ax.plot(xvals,diag[k][1:]*facs[k],clrs[k])
	ax.plot(xvals,(sum([diag[k][1:]*facs[k] for k in range(4)])*abs(kdiag))[1:],'k')
	ax.set_title('diagonal only')
	ax.set_ylim((10**-5,10**7))
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('wavevector ${nm}^{-1}$')
	ax.grid(True)		
	
	if 0: plt.tight_layout(pad=1.0,w_pad=1.0,h_pad=2.0) 
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
	plt.savefig(allsets['pickles']+'COUPLESOLN-'+str(args.passed)+'.png')
	plt.show()

#---SOLUTION LOOPS
#-------------------------------------------------------------------------------------------------------------

#---functionalized the code above which is preserved for posterity

def couplesolve(hypothesis):

	'''Function which combines the 	'hypothesize' and 'coupling_solution' sections above for looping.'''

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
	
	#---"debugging" the math here
	qs = array(meshgrid(arange(m2),arange(n2))).T
	lenscale = ms.mset.lenscale
	us = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
	qmagstd = sqrt(sum((us/array([Lx,Ly])*2*pi*lenscale)**2,axis=1))
	#---note that qdq is equal to ms.qdotq
	qdq = outer(qmagstd,qmagstd)
	if not all(qdq == ms.qdotq): raise Exception('except: problem with wavevectors')
	qd = [qdq[i,i] for i in range(len(qdq))]
	#---recap from earlier
	qmagsshift = ms.mset.lenscale*array([[sqrt(((i-m2*(i>m2/2))/((Lx)/1.)*2*pi)**2+
		((j-n2*(j>n2/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])
	xvals = reshape(qmagsshift,-1)[1:]
	#---redefine
	qdq = outer(qmagsshift,qmagsshift)
	qdqd = qdq[range(len(qdq)),range(len(qdq))]

	#---"debugging" the math here by fixing qdotq	
	termdat_test = (termlist[0]*qdq*qdq+termlist[1]*qdq+termlist[2]+termlist[3])*abs(ms.kqqp)

	diag = array([termdat_test[i,i] for i in range(len(qvals))])
	yvals = diag[1:]
	
	#---having used the diagonal of termdat_test to establish the "simple" fluctuation spectrum
	#---...we now use the blurry binner to retrieve the wavevector magnitudes after reducing dimensionality
	xvals3,yvals3,inds = blurry_binner(xvals,yvals,bin_width=bin_width)

	#---the maximum index returned by blurry_binner is equal to the number of reduced dimensions
	reddim = inds.max()
	
	#---tiny averages all magnitudes of the target data (termdat_test) with wavevectors in the i,j bin
	tiny = zeros((reddim,reddim))
	#---tinycount records the number of values that contributed to the mean
	tinycount = zeros((reddim,reddim))

	#---loop over bins and take average
	st = time.time()
	for i in range(reddim):
		status('i = '+str(i),start=st,i=i,looplen=reddim)
		for j in range(reddim):
			#---check code once
			if i == 5 and j == 3:
				#---the jjj matrix provides a matrix of combinations of positions in inds equal to i,j
				jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
				#---jj provides another way to extract the data from termdat_test via a 1D array
				jj = reshape(jjj,(product(shape(jjj)[:2]),2))
				#---we can use the 1D version and reshape it to get the "sieve" of correct values
				sieve_via_reshape = reshape(termdat_test[jj[:,0],jj[:,1]],shape(jjj)[:2])
				#---but the 1D method is redundant because we can access it directly as well
				sieve_via_direct = termdat_test[jjj[...,0],jjj[...,1]]
				#---we check that both methods agree
				reshape_passes = all(sieve_via_reshape==sieve_via_direct)
				#---and more importantly, we confirm that e.g. the last element in jjj is correct
				example_passes = all([inds[jjj[-1,-1][zz]]==[i,j][zz] for zz in range(2)])
				#---break if the test fails
				if not reshape_passes or not example_passes: raise Exception('except: code broke')				
			#---find indices i,j in the inds list and take all combinations via transposed meshgrid
			jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
			#---here we "sieve" the correct values from termdat_test
			sieve = termdat_test[jjj[...,0],jjj[...,1]]
			#---average the termdat_test values with wavevectors with indices i,j
			tiny[i,j] = mean(sieve)
			tinycount[i,j] = len(sieve)

	#---compute eigenspectrum
	fullans = linalg.eig(abs(tiny))
	prefactor = 2*product(mean(ms.mset.vecs,axis=0)[:2]/10.)
	combo = mean([fullans[1][i]*fullans[0][i] for i in range(len(fullans[1]))],axis=0)
	del ms,df
	return (xvals3[1:],prefactor*abs(combo))
	
if 'coupling_solution_loops' in routine:

	if 'loopdat' not in globals():
		loopdat = []
		for h in range(len(hypotheses)):
			status('testing hypothesis '+str(h))
			loopdat.append(couplesolve(hypotheses[h]))

	resids = array([sum(abs(loopdat[j][1][where(loopdat[j][0]<1.0)]-1.)**2)
		for j in range(len(loopdat))])
	print hypotheses[argsort(resids)[0]]
	ax = plt.subplot(111)
	alphas = array([1.0]+list(linspace(0.5,0.2,len(loopdat)-1)))[argsort(argsort(resids))]
	for l in range(len(loopdat)):
		ax.plot(loopdat[l][0],loopdat[l][1],('ro-'if alphas[l]==1. else 'k-'),
			lw=2,alpha=alphas[l],)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('wavevector ${nm}^{-1}$')
	ax.set_title('coupled modes')
	ax.grid(True)
	plt.show()
	
	

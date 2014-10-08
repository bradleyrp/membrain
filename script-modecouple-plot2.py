#!/usr/bin/python

import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---settings
fspec = 'homog' if hypothesis['kappa']['fore'] == hypothesis['kappa']['back'] else 'heterog'
combosoln_rowsum = [
	range(-10,0),
	range(1000),
	][1]
toscreen = False

if 'routine' not in globals(): routine = []

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

#---TESTING
#-------------------------------------------------------------------------------------------------------------

#---redoing termlist from scratch
if 'termdat' not in globals():

	m2,n2 = shape(ms.hqs)[1:]-0*array([1,1])
	st = time.time()
	termdat = zeros((m2*n2,m2*n2))
	for fr in range(len(ms.hqs)):
		status('fr = '+str(fr),start=st,i=fr,looplen=len(ms.hqs))
		#---note also removed trim from here
		#---note that this tests the outer function
		termdat += abs(outer(ms.hqs[fr],ms.hqs[fr]))/float(len(ms.hqs))
		if 0: termdat += outer(abs(ms.hqs[fr]),abs(ms.hqs[fr]))/float(len(ms.hqs))
		
if 1:
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
	
if 1:
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
	#---corrected method that uses the whole termdat
	tmp = array([termdat[i,i] for i in range(len(qvals))])
	yvals = tmp[1:]
	#---plot after collapsing in three different ways	
	if 1: 
		#---raw data
		print 'plotting raw data...'
		plot_spectrum1d(xvals,yvals,'',show=True,logplot=True,)	
	if 1:
		#---perfect binner
		print 'plotting data after perfect binning...'
		xvals2,yvals2 = perfect_collapser(xvals,yvals)
		plot_spectrum1d(xvals2,yvals2,'',show=True,logplot=True,)	
	if 1:
		#---blurry binner
		print 'plotting data after blurry binning...'
		xvals3,yvals3,inds = blurry_binner(xvals,yvals)
		plot_spectrum1d(xvals3,yvals3,'',show=True,logplot=True,)	
		
#---COMPUTE
#-------------------------------------------------------------------------------------------------------------

#---DIMENSIONALITY REDUCTION
if 1:
	#---working on the dimensionality reduction and checking that meshgrid works
	if 1:
		xvals3,yvals3,inds = blurry_binner(xvals,yvals)
		ttt = array(meshgrid(inds,inds)).T
		all(ttt[100,:,0]==inds[100])
		all(ttt[100,:,1]==inds)
		#---in this case we have constructed a list of all combinations of the bin indices
		jjj = array(meshgrid(where(inds==65)[0],where(inds==80)[0])).T
		jj = reshape(jjj,(product(shape(jjj)[:2]),2))
		termdat[jj[:,0],jj[:,1]]
	reddim = inds.max()
	tiny = zeros((reddim,reddim))
	tinycount = zeros((reddim,reddim))
	st = time.time()
	for i in range(reddim):
		status('i = '+str(i),start=st,i=i,looplen=reddim)
		for j in range(reddim):
			if 0:
				jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
				jj = reshape(jjj,(product(shape(jjj)[:2]),2))
				sieve = termdat[jj[:,0],jj[:,1]]
			#---I did one check, but this needs another to make sure indexing is working right
			jjj = array(meshgrid(where(inds==i)[0],where(inds==j)[0])).T
			sieve = termdat[jjj[...,0],jjj[...,1]]
			tiny[i,j] = mean(sieve)
			tinycount[i,j]  = len(sieve)
	plt.imshow(array(tiny/tinycount).T,origin='lower',interpolation='nearest',norm=mpl.colors.LogNorm())
	plt.show()


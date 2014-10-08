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

if 0: routine = ['eigenspec','plots','spec'][1:]
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
	qmagfilter = [0.001,0.5]
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

#---MAIN /// TESTING
#-------------------------------------------------------------------------------------------------------------

if 0:
	#---settings
	m2,n2 = shape(ms.mset.surf[0])
	#---show that the full_matrix is very sparse when kappa constant
	any(abs(
		sum(ms.full_matrix,axis=0)-\
		ms.full_matrix[arange(m2*n2),ms.full_matrix.argmax(axis=0)]
		)>10**-8)

if 0:		
	#---create bin indices for all positions on the matrix corresponding to unique q^2*q'^2 values
	#---gave up on doing the binning like this, but may be worth revisiting
	big_inds = (outer(ms.qmagstd**2,ms.qmagstd**2)/0.05).astype(int)
	collect = zeros((shape(ms.full_matrix)[0],shape(ms.full_matrix)[1],2))
	collect[...,0] = big_inds
	collect[...,1] = ms.full_matrix
	fmr = reshape(ms.full_matrix,-1)
	bir = reshape(big_inds,-1)

if 0:
	#---using old functions on new ms.hqs to show that fftshift on either q or v values gives almost q^-4
	hqsa = autocorr(ms.hqs,direct=1,lenscale=ms.mset.lenscale)
	t2d0 = mean((abs(array(hqsa)))**2,axis=0)
	plot_spectrum1d(reshape(ms.msc.qmagst,-1),reshape(t2d0,-1),'',show=True,logplot=True)
	plot_spectrum1d(reshape(fft.ifftshift(ms.msc.qmagst),-1),reshape(t2d0,-1),'',show=True,logplot=True)
	plot_spectrum1d(reshape(ms.msc.qmagst,-1),reshape(fft.fftshift(t2d0),-1),'',show=True,logplot=True)

if 0:
	#---using old functions on new ms.hqs with fftshift FIRST gives a much cleaner spectrum as expected
	#---note that direct = 0 gives same result as direct = 1 for some reason
	hqsa = autocorr([fft.fftshift(i) for i in ms.hqs],direct=1,lenscale=ms.mset.lenscale)
	t2d0 = mean((abs(array(hqsa)))**2,axis=0)
	plot_spectrum1d(reshape(ms.msc.qmagst,-1),reshape(t2d0,-1),'',show=True,
		logplot=True,ylims=[10**-7,10**0])
		
if 0:
	#---PROCEED AND MAKE THIS MATCH WITH NEWER "OUTER" METHOD ...
	#---using old functions on new ms.hqs with fftshift FIRST gives the spectrum with the 
	#---...ugly swoop at high q
	#---...and can be done with a re-reshaped qmagst that I calculated so now the entire 
	#---...calculation is performed
	#---...on the "new" data with the "old" functions
	hqsa = autocorr(ms.hqs,direct=1,lenscale=ms.mset.lenscale)
	t2d0 = mean((abs(array(hqsa)))**2,axis=0)
	plot_spectrum1d(reshape(reshape(ms.qmagstd,(66,66))[:-1,:-1],-1),reshape(t2d0,-1),'',
		show=True,logplot=True)

if 0:
	#---now try without the redundant high-q part
	#---...but this requires throwing out data so it's not worth it
	m2,n2 = shape(ms.mset.surf[0])
	m2,n2 = m2-1,n2-1
	hqsa = autocorr(ms.hqs,direct=1,lenscale=ms.mset.lenscale)
	t2d0 = mean((abs(array(hqsa)))**2,axis=0)
	plot_spectrum1d(reshape(reshape(ms.qmagstd,(66,66))[:-1,:-1],-1)[:m2*n2/16],reshape(t2d0,-1)[:m2*n2/16],
		'',show=True,logplot=True)

if 0:
	#---now try without the redundant high-q part
	#---...still looks spread out
	#---...re-folding the data fails probably due to the number of rows being redundant basically
	m2,n2 = shape(ms.mset.surf[0])
	m2,n2 = m2-1,n2-1
	Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]
	hqsa = autocorr(ms.hqs,direct=1,lenscale=ms.mset.lenscale)
	qmagshift = sqrt(sum((ms.usshift/(array([Lx,Ly])*ms.mset.lenscale/pi))**2,axis=1))
	t2d0 = mean((abs(array(hqsa)))**2,axis=0)
	plot_spectrum1d(reshape(reshape(qmagshift,(66,66))[:-1,:-1],-1),reshape(t2d0,-1),
		'',show=True,logplot=True)

'''
UPCOMING STEPS:
	reformulate the qmagst without the shift in order to do the autocorr without the shift
	get rid of autocorr and start matching to the hqhq from the original computation, which is on 66x66
	redo the entire process using the outer method and the nonzero values of kqqp for constant kappa
'''

def perfect_collapser(xs,ys):
	'''
	Return two arrays for a "perfect" collapse.
	'''
	xsort = array(sort(list(set(xs)))[1:])
	#inds = argmax(array([qmagstd==xsort[i] for i in range(len(xsort))]),axis=0)
	inds = argmax(array([xs==xsort[i] for i in range(len(xsort))]),axis=0)
	col = array([mean(ys[where(inds==i)]) for i in range(len(xsort))])
	return xsort,col

#---redoing termlist from scratch
if 'termdat' not in globals():

	#---original method before debugging the single-frame below
	if 0:
		#---note I removed the trim here because otherwise the non-perfect-collapsed data spread out at low q
		m2,n2 = shape(ms.hqs)[1:]-0*array([1,1])
		st = time.time()
		termdat = zeros((m2*n2,m2*n2))
		for fr in range(len(ms.hqs)):
			status('fr = '+str(fr),start=st,i=fr,looplen=len(ms.hqs))
			#---note also removed trim from here
			#---note that this tests the outer function
			termdat += abs(outer(ms.hqs[fr],ms.hqs[fr]))/float(len(ms.hqs))
			if 0: termdat += outer(abs(ms.hqs[fr]),abs(ms.hqs[fr]))/float(len(ms.hqs))
			
		#---generate the magnitudes
		Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]
		cm,cn = [int(round(i/2.-1)) for i in shape(ms.hqs)[1:]]
		qs = array(meshgrid(arange(m2),arange(n2))).T
		us = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
		qmagstd = sqrt(sum(((us-1*(us>array([cm,cn]))*\
			array([m2,n2]))/array([Lx,Ly])*2*pi*ms.mset.lenscale)**2,axis=1))
		yvals = termdat[arange(m2*n2),termdat.argmax(axis=0)]
		if 0: plot_spectrum1d(qmagstd,yvals*qmagstd**2,'',show=True,logplot=True,ylims=[10**-3,10**2])
		if 0: plot_spectrum1d(collapse_spectrum(reshape(qmagstd,(m2,n2)),reshape(qmagstd,(m2,n2))),
			collapse_spectrum(reshape(qmagstd,(m2,n2)),reshape(yvals,(m2,n2))),'',show=True,logplot=True)
		xvals2,yvals2 = perfect_collapser(qmagstd,yvals)
		#---this will plot all of the 4225 points to contrast the collapser which averages those at distinct q
		#---however note that this looks suspiciously spread out over high q so I suspect a problem
		if 0: plot_spectrum1d(reshape(qmagstd,-1),reshape(yvals,-1),'',show=True,logplot=True,ylims=[10**-3,10**2])
		#---this averages those at distinct q
		if 0: plot_spectrum1d(xvals2,yvals2,'',show=True,logplot=True,ylims=[10**-3,10**2])

	#---added back fftshift after I riffed on the old-old built-in method and saw the double swoop error
	#---this is on hold for now
	if 1:
		
		hqs = [fft.fftshift(fft.fft2(ms.mset.surf[fr][:-1,:-1])) for fr in range(len(ms.mset.surf))]
	
		m2,n2 = shape(hqs)[1:]-0*array([1,1])
		st = time.time()
		termdat = zeros((m2*n2,m2*n2))
		for fr in range(len(hqs)):
			status('fr = '+str(fr),start=st,i=fr,looplen=len(hqs))
			#---note also removed trim from here
			#---note that this tests the outer function
			termdat += abs(outer(hqs[fr],hqs[fr]))/float(len(hqs))
			if 0: termdat += outer(abs(hqs[fr]),abs(hqs[fr]))/float(len(hqs))	


	#---generate the magnitudes
	


#---? slowly debugging a single frame
if 1:

	if 0:

		#---check single frame maxima
		tmp = outer(ms.hqs[0][:-1,:-1],ms.hqs[0][:-1,:-1])
		print '(ms.msc.hqs[0]**2).max() = ...\n\t',
		print (ms.msc.hqs[0]**2).max()
		print 'outer(ms.hqs[0][:-1,:-1],ms.hqs[0][:-1,:-1]).max() = ...\n\t',
		print tmp.max()

		#---check single frame minima
		tmp = outer(ms.hqs[0][:-1,:-1],ms.hqs[0][:-1,:-1])
		print '(ms.msc.hqs[0]**2).min() = ...\n\t',
		print (ms.msc.hqs[0]**2).min()
		print 'outer(ms.hqs[0][:-1,:-1],ms.hqs[0][:-1,:-1]).min() = ...\n\t',
		print tmp.min()
	
	#---maxima and minima are equivalent, which means that the error must be downstream
	#---in this case we consider what happens to these data after they are computed
	
	#---before continuing confirm that the msc.t2d[0] scales as q4
	if 0:
		xvo,yvo = perfect_collapser(reshape(ms.msc.qmagst,-1),reshape(fft.fftshift(ms.msc.t2d[0]),-1))
		plot_spectrum1d(xvo,yvo,'',show=True,logplot=True,)

		#---slope was only -3.2 so now checking without perfect collapser	
		xvo2,yvo2 = (reshape(fft.fftshift(ms.msc.qmagst),-1),reshape(ms.msc.t2d[0],-1))
		plot_spectrum1d(xvo2,yvo2,'',show=True,logplot=True,)
		#---this shows a bit of spread which might be why it's 3.2 and not 4
		#---it also requires the fftshift above or you get the spike in the wrong place
	
		#---manually calculating the spectrum from the old ms.msc.hqs data
		yvo3 = reshape(abs(mean(ms.msc.hqs[:,:,:]**2,axis=0)),-1)[1:]
		qmags = ms.mset.lenscale*array([[sqrt(((i)/((Lx)/1.)*2*pi)**2+((j)/((Ly)/1.)*2*pi)**2)
			for j in range(0,n2)] for i in range(0,m2)])
		xvo3 = reshape(qmags,-1)[1:]
		plot_spectrum1d(xvo3,yvo3,'',show=True,logplot=True,)

	if 0:
	
		#---using the built-in undulation calculator instead of the one from ms.msc
		#---note that this is the oldest and most reliable way to calculate
		#---? does this indicate errors in ms.msc?
		ys = reshape(ms.mset.undulate_hqhq2d,-1)
		xs = reshape(ms.mset.undulate_qmag2d,-1)
		plot_spectrum1d(xs,ys,'',show=True,logplot=True,)
		#---this is very close to 4 and doesn't show any weird hi-q behavior
		#---hence this will be my reference for the rest of this exercise

	'''
	steps
		recap the built-in undulation calculator
			compute fft with trimmed redundant rows and send to undulate_raw
			undulate_hqhq2d = mean(array(1.*(abs(array(undulate_raw))/double(m*n))**2),axis=0)
			spec_center = [int(i/2) for i in shape(self.undulate_raw)[1:]]
			qs = [(sqrt((2*pi*(array(xv)-spec_center[0])/lxy[f][0])**2+
				(2*pi*(array(yv)-spec_center[1])/lxy[f][1])**2)) for f in range(len(lxy))]
			undulate_qmag2d = mean(qs,axis=0)
		compare outer vs reliable method
	'''	

	if 0:	

		#---remake of the built-in undulation calculator
		m,n = m2-1,n2-1
		lxy = ms.mset.vecs
		yv,xv = meshgrid(range(n),range(m))
		undulate_raw = [fft.fftshift(fft.fft2(ms.mset.surf[fr][:-1,:-1])) for fr in range(len(ms.mset.surf))]
		undulate_hqhq2d = mean(array(1.*(abs(array(undulate_raw))/double(m*n))**2),axis=0)
		spec_center = [int(i/2) for i in shape(undulate_raw)[1:]]
		qs = ms.mset.lenscale*array([(sqrt((2*pi*(array(xv)-spec_center[0])/lxy[f][0])**2+
			(2*pi*(array(yv)-spec_center[1])/lxy[f][1])**2)) for f in range(len(lxy))])
		undulate_qmag2d = mean(qs,axis=0)
		ys = reshape(undulate_hqhq2d,-1)[1:]
		xs = reshape(undulate_qmag2d,-1)[1:]
		plot_spectrum1d(xs,ys,'',show=True,logplot=True,)	
		#---this perfectly captures the original code


		#---riff with no fftshift
		m,n = m2-1,n2-1
		lxy = ms.mset.vecs
		yv,xv = meshgrid(range(n),range(m))
		undulate_raw = [fft.fft2(ms.mset.surf[fr][:-1,:-1]) for fr in range(len(ms.mset.surf))]
		undulate_hqhq2d = mean(array(1.*(abs(array(undulate_raw))/double(m*n))**2),axis=0)
		spec_center = [int(i/2) for i in shape(undulate_raw)[1:]]
		qs = ms.mset.lenscale*array([(sqrt((2*pi*(array(xv))/lxy[f][0])**2+
			(2*pi*(array(yv))/lxy[f][1])**2)) for f in range(len(lxy))])
		undulate_qmag2d = mean(qs,axis=0)
		ys = reshape(undulate_hqhq2d,-1)[1:]
		xs = reshape(undulate_qmag2d,-1)[1:]
		plot_spectrum1d(xs,ys,'',show=True,logplot=True,)
		#---! this gives you the large double-near swoop at the end
		
		#---riff with no fftshift, no trimming
		m,n = m2,n2
		lxy = ms.mset.vecs
		yv,xv = meshgrid(range(n),range(m))
		undulate_raw = [fft.fft2(ms.mset.surf[fr][:,:]) for fr in range(len(ms.mset.surf))]
		undulate_hqhq2d = mean(array(1.*(abs(array(undulate_raw))/double(m*n))**2),axis=0)
		spec_center = [int(i/2) for i in shape(undulate_raw)[1:]]
		qs = ms.mset.lenscale*array([(sqrt((2*pi*(array(xv))/lxy[f][0])**2+
			(2*pi*(array(yv))/lxy[f][1])**2)) for f in range(len(lxy))])
		undulate_qmag2d = mean(qs,axis=0)
		ys = reshape(undulate_hqhq2d,-1)[1:]
		xs = reshape(undulate_qmag2d,-1)[1:]
		plot_spectrum1d(xs,ys,'',show=True,logplot=True,)
		#---! same result with almost no change at all
	

#---other single-frame testing notes
if 0:

	#---? check magnitude maxima
	hqsa = autocorr(ms.msc.hqs[:1],direct=1,lenscale=ms.mset.lenscale)
	t2d0frame = mean((abs(array(hqsa)))**2,axis=0)



	'''
	>>> t2d0frame.max()
	0.04117415113807435
	>>> abs(tmp.max())
	42.0701368717043
	>>>  abs(tmp).max()
	42.0701368717043
	>>> abs((ms.msc.hqs[0]**2).max())
	42.070136871704285
	'''
	'''
	>>> abs(tmp.min())
	36.9976291261368
	>>> abs((ms.msc.hqs[0]**2).min())
	36.9976291261368
	>>> abs(tmp).min()
	2.0779449278478598e-26
	>>> abs(ms.msc.hqs[0]**2).min()
	2.0779449278478598e-26
	'''
	'''
	>>> sort(reshape(abs(ms.msc.hqs[0]**2),-1))[1]
	1.3623469053561301e-08
	>>> sort(reshape(abs(tmp),-1))[1]
	1.6825224640265576e-17
	'''

	
	
	
			
#---DEVELOPMENT
#-------------------------------------------------------------------------------------------------------------

if hasattr(ms,'fullans'): 
	
	#---first attempt to manipulate the diagonalized matrix

	#---combine some rows to look at the coupling
	fullans = ms.fullans
	n2,m2 = ms.m2,ms.n2
	#---average some solutions by eigenvalue magnitude
	rowcombo = mean([array([[fullans[1][fr][j+(i-1)*m2] for j in range(n2)] 
		for i in range(m2)])*fullans[0][fr]**2 for fr in combosoln_rowsum],axis=0)
	raws = [fullans[0][i]*reshape(fullans[1][i],(m2,n2)) for i in combosoln_rowsum]
	combosoln = mean(raws,axis=0)

if 'plots_OLD' in routine:

	#---loop over plots
	for logplot in [True,False][:1]:
		for data in [
			[ms.kqqp,'kqqp'],
			[ms.kfield,'kfield'],
			[ms.full_matrix,'full_matrix'],
			[ms.qdotq**2*termlist[0],'qqhqhq'],
			[termlist[0],'hqhq'],
			[termlist[1],'hqc0qp'],
			][:]:
			print 'data = '+data[1]+' logplot = '+str(logplot)
			plot_spectrum(data[0],callsign+'-'+fspec+'-term-'+data[1]+'-',logplot=logplot,show=toscreen)

if 'plots_dev' in routine:

	#---DEPRECATED

	#---plot eigenspectrum
	plot_spectrum1d(
		range(len(fullans[0]))[:100],
		abs(fullans[0])[:100],
		callsign+'-'+fspec+'EIGENSPEC',show=toscreen)

	#---SPEC
	plot_spectrum1d(
		ms.qmagunshift[1:],
		sum(abs(ms.full_matrix),axis=1)[1:],
		callsign+'-'+fspec+'SPEC',show=toscreen)

	#---SPEC2
	plot_spectrum1d(
		qmagunshift[1:],
		sum(abs(fmuns),axis=1)[1:],
		callsign+'-'+fspec+'SPEC2',show=toscreen)
	
	#---SPECTEST
	plot_spectrum1d(
		ms.qmagunshift[1:],
		(sum(termlist[0],axis=1)*ms.qmagunshift**2)[1:],
		callsign+'-'+fspec+'SPECTEST',show=toscreen)

	#---SPECTESTb
	plot_spectrum1d(
		reshape(fft.ifftshift(reshape(ms.qmagunshift,(66,66))),-1)[1:],
		(sum(termlist[0],axis=1)*ms.qmagunshift**2)[1:],
		callsign+'-'+fspec+'SPECTESTb',show=toscreen)

	#---SPECTEST2
	plot_spectrum1d(
		ms.qmagunshift[1:],
		sum(abs(ms.full_matrix),axis=1)[1:],
		callsign+'-'+fspec+'SPECTEST2',show=toscreen)

	#---SPECTEST3
	qmagshift = sqrt(sum((ms.usshift/(array([Lx,Ly])*ms.mset.lenscale/pi))**2,axis=1))
	plot_spectrum1d(
		sum(abs(termlist[0]),axis=1)[1:],
		reshape(ms.msc.qmagst,-1)[1:],
		callsign+'-'+fspec+'SPECTEST3',show=toscreen)

	#---SPEC3
	plot_spectrum1d(
		ms.qmagunshift[1:],
		(sum(abs(ms.full_matrix),axis=1)*ms.qmagunshift**2)[1:],
		callsign+'-'+fspec+'SPEC3',show=toscreen)

	#---spec4undone
	complq = sum(array(meshgrid(arange(m2),1j*arange(n2))).T,axis=-1)
	plot_spectrum1d(
		ms.qmagunshift[1:],
		(sum(abs(ms.full_matrix),axis=1)*ms.qmagunshift**2)[1:],
		callsign+'-'+fspec+'SPEC4',show=toscreen)

if 'demo_indices' in routine:

	print 'demonstrating outer product indexing by doing the lookup three different ways'
	dddn,dddm = 3,3
	a = array([[2,4,6],[8,10,12],[14,16,18]])
	b = array([[2,4,6],[8,10,12],[14,16,18]])+1
	i,j = 3,6
	print outer(a,b)
	print outer(a,b)[i,j]
	print reshape(a,-1)[i]*reshape(b,-1)[j]
	print a[i/dddm,i%dddm]*b[j/dddn,j%dddn]

if 'further_development' in routine:

	#---single frame hqhq test
	#---INCLUDES REDUNDANT EDGE, NOT FFTSHIFTED
	if 1:

		trim = False
		shift = False

		m2,n2 = shape(ms.mset.surf[0])
		if trim: m2,n2 = m2-1,n2-1
		cm,cn = [int(round(i/2.-1)) for i in shape(ms.mset.surf[0])]
		Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]

		frames_to_sample = 100
		st = time.time()
		termdat = zeros((m2*n2,m2*n2))
		hqs = [None for i in range(frames_to_sample)]
		for fr in range(frames_to_sample):
			status('fr = '+str(fr),start=st,i=fr,looplen=frames_to_sample)
			hqs[fr] = fft.fft2(ms.mset.surf[fr][:(-1 if trim else None),:(-1 if trim else None)])
			if shift: hqs[fr] = fft.fftshift(hqs[fr])
			termdat += abs(outer(hqs[fr],hqs[fr]))
		termdat = termdat / float(frames_to_sample)

		hqhq1abs = termdat

		if not shift:
			qmags = ms.mset.lenscale*array([ (i)/((Lx)/1.)*2*pi+1j*(j)/((Ly)/1.)*2*pi 
				for j in range(0,n2) for i in range(0,m2)])
			qmags_opp = ms.mset.lenscale*array([ (j)/((Lx)/1.)*2*pi+1j*(i)/((Ly)/1.)*2*pi 
				for j in range(0,n2) for i in range(0,m2)])
		else:
			qmags = ms.mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi 
				for j in range(0,n2) for i in range(0,m2)])
			qmags_opp = ms.mset.lenscale*array([ (j-cn)/((Lx)/1.)*2*pi+1j*(i-cm)/((Ly)/1.)*2*pi 
				for j in range(0,n2) for i in range(0,m2)])

		q4 = outer(abs(qmags**2),abs(qmags**2))
		qmagsq = sqrt(abs(outer(qmags,qmags)))

	#---MOAR plots
	if 0: 
		plot_spectrum(hqhq1abs,'',show=True,logplot=True)
		plot_spectrum1d(
			qmags,
			hqhq1abs[arange(m2*n2),arange(m2*n2)],
			'',show=True,logplot=True)
		plot_spectrum1d(
			reshape(qmagsq,-1),
			reshape(hqhq1abs,-1),
			'',show=True,logplot=True,ylims=[10**1,10**10])
		plot_spectrum(termdat,'',show=True,logplot=True)
		plot_spectrum(qmagsq,'',show=True,logplot=True)

	plot_spectrum1d(
		reshape(qmagsq,-1),
		reshape(abs(termlist[0]),-1),
		callsign+'-'+fspec+'SPECTEST10',show=True)

	plot_spectrum1d(
		reshape(qmagsq,-1),
		reshape(abs(ms.full_matrix),-1),
		callsign+'-'+fspec+'SPECTEST11',show=False)
	
	#---SPECTEST2
	plot_spectrum1d(
		ms.qmagunshift[1:],
		sum(abs(ms.full_matrix),axis=1)[1:],
		callsign+'-'+fspec+'SPECTEST2',show=toscreen)
	
	#---SPECTEST2x
	plot_spectrum1d(
		qmagsq[arange(m2*n2),arange(m2*n2)],
		sum(abs(ms.full_matrix),axis=1),
		callsign+'-'+fspec+'SPECTEST2x',show=toscreen)

	#---SPECTEST2y
	plot_spectrum1d(
		qmagsq[arange(m2*n2),arange(m2*n2)],
		abs(ms.full_matrix)[arange(m2*n2),arange(m2*n2)],
		callsign+'-'+fspec+'SPECTEST2y',show=toscreen)
	
	#---SPECTEST3h
	qeqnegq = array(where(abs(ms.kqqp)>10**-8))
	plot_spectrum1d(
		qmagsq[qeqnegq[0],qeqnegq[1]],
		abs(ms.full_matrix)[qeqnegq[0],qeqnegq[1]],
		callsign+'-'+fspec+'SPECTEST3h',show=toscreen)
	
#---NOTES
#-------------------------------------------------------------------------------------------------------------

'''
NOTES 3

check that the homogeneous full_matrix has argmax along one dimension equal to the argmax of kqqp:
	all(argmax(ms.kqqp,axis=0)==argmax(ms.full_matrix,axis=0))
	
don't shift anything
	qmagunshift = sqrt(sum((self.us/(array([Lx,Ly])*self.mset.lenscale/pi))**2,axis=1))
	qdotqun = outer(qmagunshift,qmagunshift)
	hqs_hqps=termlist[0]
	hqs_c0qps=termlist[1]
	c0qs_hqps=termlist[2]
	c0qs_c0qps=termlist[3]
	fmuns = (qdotqun*qdotqun*hqs_hqps-qdotqun*hqs_c0qps-qdotqun*c0qs_hqps+c0qs_c0qps)*real(self.kqqp)
	plot_spectrum(fmuns.T,'',show=True,logplot=True)
	
confirm that my new, unshifted qdotq is really the correct magnitude
	>>> array(where(sort(reshape(ms.qdotq,-1))!=0)).T[0]
	array([8711])
	sqrt(sort(reshape(ms.qdotq,-1))[8711])
	>>> sqrt(sort(reshape(ms.qdotq,-1))[8711])
	>>> sort(abs(ms.qmags_nc))[1]
	0.095811928046317107
	
working on comparing new ModeSum method to old ModeCouple and autocorr code to troubleshoot spectrum
	first I added the old code for accessing the ModeCouple module (deleted a few days ago) to compare
	hqsa = autocorr(ms.msc.hqs,direct=0,lenscale=ms.mset.lenscale)
	
	tracing things back from the script-coupling-adv-batch.py script
		yvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.t2d[0]))
		xvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.qmagst))
		ax=plt.subplot(111);ax.plot(xvals1dspecold,yvals1dspecold);ax.set_xscale('log');ax.set_yscale('log');plt.show()
		where does t2d come from?
			mean((abs(array(hqsa)))**2,axis=0)
		where does hqsa come from?
			hqsa = autocorr(ms.msc.hqs,direct=0,lenscale=ms.mset.lenscale)
		can I do this outside ModeCouple?
			hqsa = autocorr(ms.msc.hqs,direct=0,lenscale=ms.mset.lenscale)				
			t2drecap = collapse_spectrum(ms.msc.qmagst,mean((abs(array(hqsa)))**2,axis=0))
		confirmed with 
			all(t2drecap == ms.msc.t2d[0])
		using the current hqs from ms and not ms.msc
			hqsanew = autocorr(ms.hqs,direct=0,lenscale=ms.mset.lenscale)
			t2drecapnew = array(collapse_spectrum(ms.msc.qmagst,mean((abs(array(hqsanew)))**2,axis=0)))
		does this have the discrepancy?
			yes:
				>>> t2drecap.max()/t2drecapnew.max()
				48.433103481381302
		is this due to averaging that happens in collapse_spectrum?
			mean((abs(array(hqsa)))**2,axis=0).max()/mean((abs(array(hqsanew)))**2,axis=0).max()
		or is it due to the hqsa step?
			>>> abs(array(ms.msc.hqs[0]).max()),abs(array(ms.hqs[0]).max())
			(2.0403309681286368, 2.0291414721027796)
			>>> abs(array(ms.msc.hqs[0]).max()),abs(array(ms.hqs[0][:-1,:-1]).max())
			(2.0403309681286368, 1.3087569206388452)
			>>> shape(ms.hqs[0][:-1,:-1])
			(65, 65)
			>>> shape(ms.msc.hqs[0])
			(65, 65)
	determine what makes the good spectrum
		supposed steps for the good spectrum
			1. hqs[fr] = fftwrap(mset.surf[fr])/double(m*n)
				equivalent to fft.fftshift(fft.fft2(array(ms.mset.surf[fr])[:-1,:-1]))
				note that
					>>> fft.fftshift(fft.fft2(array(ms.mset.surf[fr])[:-1,:-1])).max()
					(8283.2086092826175+2387.4091985462255j)
			2. hqsa = autocorr(hqs,direct=0,lenscale=ms.mset.lenscale)
				executed 
					hqexnew = (ms.fftwrap(array(ms.mset.surf[fr]))/double(product(ms.mset.griddims)))
					hqexold = array([fft.fftshift(fft.fft2(array(ms.mset.surf[fr])[:-1,:-1])) for fr in range(500)])
					hqsaold = autocorr(hqexold,direct=0,lenscale=ms.mset.lenscale)	
					t2d0old = mean((abs(array(hqsaold)))**2,axis=0)
					plt.imshow(t2d0old.T,interpolation='nearest',origin='lower');plt.show()
				which looks like a normal spectrum
				comparing the hqs and hqsa examples from the "old" style FFT convention shows 
					a factor of 10 difference which must be due to the autocorr function
					>>> abs(hqexold[0].max())
					8620.3983403434904
					>>> abs(hqsaold[0].max())
					862.03983403434916
		what I'm doing on the new spectrum
			1. hqs[fr] = ms.fftwrap(array(ms.mset.surf[fr]))/double(product(ms.mset.griddims))
				equivalent to 
					fft.fft2(array(ms.mset.surf[fr]))
				note that
					>>> (ms.fftwrap(array(ms.mset.surf[fr]))/double(product(ms.mset.griddims))).max()
					(1.95239482789925+0.55278345651773309j)
					>>> fft.fft2(array(ms.mset.surf[fr])).max()
					(8504.6318703291327+2407.9247365912452j)
					>>> fft.fft2(array(ms.mset.surf[fr][:-1,:-1])).max()
					(8283.2086092826175+2387.4091985462255j)
			2. hqsa ... ?
		so step 1 where we do the fft on the data with the final row and column dropped 
			gives the EXACT same maximum value of the outgoing data 
			whereas the version I'm using is slightly higher
		another try
			>>> mean((abs(array(autocorr(ms.msc.hqs,direct=0,lenscale=ms.mset.lenscale))))**2,axis=0).max()
			0.23960450528727523
			>>> mean((abs(array(ms.msc.hqs)))**2,axis=0).max()23.960450528727538
			23.960450528727538
			>>> mean((abs(array(ms.msc.hqs)/ms.mset.lenscale))**2,axis=0).max()
			0.23960450528727523
		try a short computation to check on things
			b0 = array([outer(ms.hqs[fr],ms.hqs[fr]) for fr in range(10)])
			b1 = mean(b0,axis=0)
			c0 = array([fft.fftshift(fft.fft2(array(ms.mset.surf[fr])[:-1,:-1])) for fr in range(10)])
			c1a = autocorr(c0,direct=0,lenscale=ms.mset.lenscale)	
			c1 = mean((abs(array(c1a)))**2,axis=0)
		RETHINK BEFORE COFFEE
		REBOOT LATER IN THE EVENING
		start with a nice block of code that recaps the hqhq really well
			yvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.t2d[0]))
			xvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.qmagst))
			ax=plt.subplot(111);ax.plot(xvals1dspecold,yvals1dspecold);ax.set_xscale('log');ax.set_yscale('log');plt.show()
		new plan: find the maximum of this, find the q values it corresponds to, and make sure the termlist[0] outerproduct has the same value in the same position
			plot the t2d[0] value:
				plot_spectrum(ms.msc.t2d[0].T,'',show=True,logplot=True)
				but this is done via fftshift so it's hard to really tell what's going on
			reconstitute the t2d without the fftshift
				USE THE FOLLOWING CODE AND MODIFY THE FLAGS IN ModeCouple
					yvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.t2d[0]))
					xvals1dspecold = array(collapse_spectrum(ms.msc.qmagst,ms.msc.qmagst))
					ax=plt.subplot(111);ax.plot(xvals1dspecold,yvals1dspecold);ax.set_xscale('log');ax.set_yscale('log');plt.show()
					plot_spectrum(ms.msc.qmagst,'',logplot=True,show=True)
				with redundant = True, shift = False
					...gives qmagst relative to the center point
					...and a FLAT spectrum t2d[0] which looks pretty sweet
				with redundant = True, shift = True
					...gives what looks like a nice q**-4 scaling
					...gives qmagst relative to center point
				with redundant = False, shift = True
					...steep 1d spectrum !
					...same qmagst
				with both False
					...flat
				so looking only at msc, if you don't do the shift, then things look flat, which they shouldn't befspec
					
'''			

'''
NOTES 1

m3,n3=m2-1,n2-1
inds = array(where(ms.kqqp>1)).T[1]
[ms.qs[(i)/n3,(i)%n3] for i in inds]
ms.usshift[inds[0]],ms.usshift[inds[1]]
[ms.qs[i/n2,i%n2] for i in j for j in [[1386,2970][1539,2883]]]
sqrt(sum((ms.usshift/(array([Lx,Ly])*ms.mset.lenscale))**2,axis=1))
unravel_index(termlist[0].argmax(),shape(termlist[0]))
#---useful to confirm right distance metric AND right reshape behavior
plot_spectrum(reshape(qmagshift,(m2,n2)).T,'',show=True,logplot=False)
plot_spectrum(outer(qmagshift,qmagshift).T,'',show=True,logplot=False)
'''

'''
NOTES 2

plot_spectrum(real((ms.qmags*ms.qmags)*(ms.qmags*ms.qmags)*termlist[0]).T,'',show=True,logplot=True)
plot_spectrum(real(termlist[0]).T,'',show=True,logplot=True)
plot_spectrum(real(termlist[0]).T,'',show=True,logplot=False)
plot_spectrum(ms.kqs.T,'',show=True,logplot=False)
plot_spectrum(fft.fftshift(termlist[0]).T,'',show=True,logplot=True)
plot_spectrum(fft.fftshift(termlist[1]).T,'',show=True,logplot=True)
'''

'''

####################################
The final copy of script-modecouple-plot.py that helped me match the original q4 scaling to the
new outer method. Recapped here for posterity
####################################

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
	#Temporary spectrum plotting program.
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
	#Temporary spectrum plotting program.
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
	qmagfilter = [0.001,0.5]
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

#---MAIN /// TESTING
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
		
if 0:
	#---remake of the built-in undulation calculator
	#---...which perfectly captures the original code and correct scaling
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
	plot_spectrum1d(xs,ys,'',show=True,logplot=True,)	
	

	#---new version from outer to confirm that outer works to compute the fluctuations
	Lx,Ly = mean(ms.mset.vecs,axis=0)[0:2]		
	qmags = ms.mset.lenscale*array([[sqrt(((i)/((Lx)/1.)*2*pi)**2+((j)/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])	
	qvals = reshape(qmags,-1)


#---a test
if 1:
	u,v = 200,200
	print outer(ms.hqs[0],ms.hqs[0])[u][v]
	print ms.hqs[0][u/n2][u%n2]**2
	print abs(mean([abs(ms.hqs[fr][u/n2][u%n2])**2 for fr in range(len(ms.mset.surf))]))
	print termdat[u,u]
	
#---a test, after setting trim = 0 on the remake above
#---...and since using trim made this one fail to find the 4.24e-26 value in the undulate_hqhq2d
#---...we conclude that trimming drops that very-close-to-zero value
if 1:
	tmp = array([termdat[i,i] for i in range(len(qvals))])
	for i in [0,1]:
		print sort(reshape(undulate_hqhq2d,-1))[i]
		print sort(tmp)[i]

#---visualize the right scaling but with the swoop because I haven't PBC shifted the really big wavevectors
if 0: plot_spectrum1d(qvals[1:],tmp[1:],'',show=True,logplot=True,)

#---struggling to do the fold-over properly here, so I'm tweaking my calculation of qmagsshift and then 
#---...viewing it with imshow to make sure id gives that plus-sign that indicates the correct distances
if 0:
	qmagsshift = ms.mset.lenscale*array([[sqrt(((i-m2*(i>m2/2))/((Lx)/1.)*2*pi)**2+
		((j-n2*(j>n2/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n2)] for i in range(0,m2)])
	plot_spectrum1d(reshape(qmagsshift,-1)[1:],tmp[1:],'',show=True,logplot=True,)
	plot_spectrum1d(reshape(fft.fftshift(qmags),-1)[1:],
		reshape(fft.fftshift(reshape(tmp,(m2,n2))),-1)[1:],'',show=True,logplot=True,)
'''

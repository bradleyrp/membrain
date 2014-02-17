#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

import scipy.ndimage

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---settings
xyzform = 'rect'
nbase = 22
length = None
qmagfilter=[10**-10,10**6]

#---plan
analysis_descriptors = [
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-anis-22-run1-rep-0-1500-2000',True,'',False),
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-bare-22-run1-rep-0-1500-2000',True,'',True),
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',False),
	('/store-delta/compbio/mesoscale-v2002/t3-anis-22-c0-0.05/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t3-anis-22-run1-rep-0-1500-2000',True,'',False),
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t4-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-bare-22-run1-rep-0-1500-2000',True,'',True)]
ad = analysis_descriptors[3]
vtudir,start,end,testname,ismeso,pklname,isbare = ad
barecompare = True
baresys = 4
removeavg = True
fitlims = [4,2]
forcekappa = True
mdcompare = True
mdpkl = 'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl'

#---method 
'''	   |---------load
       ||--------calc
       |||-------plot terms
       ||||------plot compare hqhq
       |||||-----plot compare hqhqq4
       ||||||----plot spectrum, 1D
       |||||||---plot phase angles '''
seq = '0000000'

#---settings
showplots = True
sskip = 4 #---sets the xy ticks spacing  in units of lenscale on any 2D spectra
clrs = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,3,7,2,1,4]]
cmap = mpl.cm.jet

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def fftwrap(dat,redundant=1):
	'''This function wraps the standard discrete FFT for a system with possible-redundant rows.'''
	trim = -1 if redundant == 1 else None 
	return fft.fftshift(fft.fft2(array(dat)[:trim,:trim]))
	
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

def calculate_mode_coupling(ms):
	'''Moved the entire mode-coupling calculation to a single function so it can be repeated easily.'''
	ms.lenscale = lenscale
	ms.calculate_undulations(removeavg=removeavg,fitlims=fitlims,forcekappa=forcekappa)
	grid = ms.griddims
	m,n = grid[0]-1,grid[1]-1
	hqs = [fftwrap(ms.surf[i])/double(m*n) for i in range(len(ms.surf))]
	cqs = [fftwrap(lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
	Lx,Ly = mean(ms.vecs,axis=0)[0:2]
	cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
	qmags = lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
		for j in range(0,n)] for i in range(0,m)])
	hqsa = autocorr(hqs,direct=0,lenscale=lenscale)
	hqsb = autocorr(hqs,direct=1,lenscale=lenscale)
	center = [cm,cn]
	cqsa = autocorr(cqs,direct=0,lenscale=lenscale)
	cqsb = autocorr(cqs,direct=1,lenscale=lenscale)
	qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
	qmagst_err = np.std(qmagst,axis=0)
	mt,nt = shape(hqsa[0])
	t0 = mean((abs(array(hqsa)))**2,axis=0)
	t1 = mean(abs(hqsa)*abs(cqsb),axis=0)
	t2 = mean(abs(cqsa*hqsb),axis=0)
	t3 = mean(abs(cqsa*cqsb),axis=0)
	t0e = std((abs(array(hqsa)))**2,axis=0)
	t1e = std(abs(hqsa)*abs(cqsb),axis=0)
	t2e = std(abs(cqsa*hqsb),axis=0)
	t3e = std(abs(cqsa*cqsb),axis=0)
	t0spec = array([[qmagst[i,j],t0[i,j]] for j in range(nt) for i in range(mt)])
	t1spec = array([[qmagst[i,j],t1[i,j]] for j in range(nt) for i in range(mt)])
	t2spec = array([[qmagst[i,j],t2[i,j]] for j in range(nt) for i in range(mt)])
	t3spec = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])
	area = double(mean([ms.vec(i)[0]*ms.vec(i)[1] for i in ms.surf_index])/ms.lenscale**2)
	scalefac = ms.undulate_kappa*area
	termsum = t0*qmagst**4-t1*qmagst**2-t2*qmagst**2+t3
	extraterms = -1*t1*qmagst**2-t2*qmagst**2+t3
	termsumspec = array([[qmagst[i,j],scalefac*termsum[i,j]] for j in range(nt) for i in range(mt)])
	termsumspec_err = array([[qmagst[i,j],
		scalefac*sqrt((t0e[i,j]*qmagst[i,j]**4)**2+(2*qmagst[i,j]**2*t1e[i,j])**2+(t3e[i,j])**2)] 
		for j in range(nt) for i in range(mt)])
	termsumspec_err[termsumspec_err>=termsumspec] = termsumspec[termsumspec_err>=termsumspec]*.999999
	t0spec2 = array([[qmagst[i,j],t0[i,j]*qmagst[i,j]**4] for j in range(nt) for i in range(mt)])
	t1spec2 = array([[qmagst[i,j],t1[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t2spec2 = array([[qmagst[i,j],t2[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t3spec2 = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
if int(seq[0]) or mset.surf == []:
	lenscale = 1.0
	if mdcompare:
		print 'loading md comparison'
		msetmd = unpickle(pickles+mdpkl)
		#---compute hypothetical curvature field for the MD simulation
		vecs = mean(mset.vecs,axis=0)
		m,n = msetmd.griddims
		getgrid = array([[[i,j] for j in linspace(0,vecs[1],n)] for i in linspace(0,vecs[0],m)])
		params = [0,0.05,vecs[0]/2.,vecs[1]/2.,5,5,0]
		c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1])
		    for j in range(n)] for i in range(m)])
	if ismeso:
		c0sraw = array(mset.load_points_vtu(vtudir,extra_props='induced_cur',
			start=start,end=end,nbase=nbase,lenscale=lenscale))[:,0]
		mset.surfacer()
		c0s = mset.surfacer_general(c0sraw)
		#---note that you have to fix the scaling here. lenscale is used in rezipgrid, problem for c0
	else:
		mset = unpickle(pickles+pklname)
	if barecompare:
		vtudirb,startb,endb,testnameb,ismesob,pklnameb,isbareb = analysis_descriptors[baresys]
		mset2 = MembraneSet()
		mset2.load_points_vtu(vtudirb,start=startb,end=endb,nbase=nbase,lenscale=lenscale)
		mset2.surfacer()
	#---infer the a0 factor to match system sizes
	lenscale = max(mean(mset2.vecs,axis=0))/(max(mean(msetmd.vecs,axis=0))/msetmd.lenscale)
	print 'reset lenscale = '+str(lenscale)
	mset.lenscale = lenscale
	if barecompare:
		mset2.lenscale = lenscale
		
#---calculate mode couplings
if int(seq[1]):
	mset.lenscale = lenscale
	mset.calculate_undulations(removeavg=removeavg,fitlims=fitlims,forcekappa=forcekappa)
	grid = mset.griddims
	m,n = grid[0]-1,grid[1]-1
	hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
	cqs = [fftwrap(lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
	Lx,Ly = mean(mset.vecs,axis=0)[0:2]
	cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
	qmags = lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
		for j in range(0,n)] for i in range(0,m)])
	hqsa = autocorr(hqs,direct=0,lenscale=lenscale)
	hqsb = autocorr(hqs,direct=1,lenscale=lenscale)
	center = [cm,cn]
	cqsa = autocorr(cqs,direct=0,lenscale=lenscale)
	cqsb = autocorr(cqs,direct=1,lenscale=lenscale)
	qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
	qmagst_err = np.std(qmagst,axis=0)
	mt,nt = shape(hqsa[0])
	t0 = mean((abs(array(hqsa)))**2,axis=0)
	t1 = mean(abs(hqsa)*abs(cqsb),axis=0)
	t2 = mean(abs(cqsa*hqsb),axis=0)
	t3 = mean(abs(cqsa*cqsb),axis=0)
	t0e = std((abs(array(hqsa)))**2,axis=0)
	t1e = std(abs(hqsa)*abs(cqsb),axis=0)
	t2e = std(abs(cqsa*hqsb),axis=0)
	t3e = std(abs(cqsa*cqsb),axis=0)
	t0spec = array([[qmagst[i,j],t0[i,j]] for j in range(nt) for i in range(mt)])
	t1spec = array([[qmagst[i,j],t1[i,j]] for j in range(nt) for i in range(mt)])
	t2spec = array([[qmagst[i,j],t2[i,j]] for j in range(nt) for i in range(mt)])
	t3spec = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])
	area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
	scalefac = mset.undulate_kappa*area
	termsum = t0*qmagst**4-t1*qmagst**2-t2*qmagst**2+t3
	extraterms = -1*t1*qmagst**2-t2*qmagst**2+t3
	termsumspec = array([[qmagst[i,j],scalefac*termsum[i,j]] for j in range(nt) for i in range(mt)])
	termsumspec_err = array([[qmagst[i,j],
		scalefac*sqrt((t0e[i,j]*qmagst[i,j]**4)**2+(2*qmagst[i,j]**2*t1e[i,j])**2+(t3e[i,j])**2)] 
		for j in range(nt) for i in range(mt)])
	termsumspec_err[termsumspec_err>=termsumspec] = termsumspec[termsumspec_err>=termsumspec]*.999999
	t0spec2 = array([[qmagst[i,j],t0[i,j]*qmagst[i,j]**4] for j in range(nt) for i in range(mt)])
	t1spec2 = array([[qmagst[i,j],t1[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t2spec2 = array([[qmagst[i,j],t2[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t3spec2 = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])
	if barecompare:
		mset2.lenscale = lenscale
		mset2.calculate_undulations(removeavg=removeavg,fitlims=fitlims,forcekappa=forcekappa)
		area = double(mean([mset2.vec(i)[0]*mset2.vec(i)[1] for i in mset2.surf_index])/mset2.lenscale**2)
		scalefacbare = mset2.undulate_kappa*area
		specbare = mset2.undulate_spec1d
		energy_bare = array([[specbare[j,0],specbare[j,0]**4*specbare[j,1]*scalefacbare] 
			for j in range(len(specbare))])
		t0_bare = autocorr([mset2.undulate_hqhq2d])[0]
	if mdcompare:
		area = double(mean([msetmd.vec(i)[0]*msetmd.vec(i)[1] for i in msetmd.surf_index])/msetmd.lenscale**2)
		scalefac = msetmd.undulate_kappa*area
		specmd = msetmd.undulate_spec1d
		energy_md = array([[specmd[j,0],specmd[j,0]**4*specmd[j,1]*scalefac] 
			for j in range(len(specmd))])
		calculate_mode_coupling(msetmd)

	
#---view individual contributions to the energy terms in two dimensions
if int(seq[2]):
	islognorm = True
	vmin = 10**-10
	vmax = 10**0
	fig = plt.figure(figsize=(8,8))
	gs = gridspec.GridSpec(2,2)
	axes = []
	ax = plt.subplot(gs[0,0])
	axes.append(ax)
	ax.set_title('term 1')
	ax.imshow(array(t0).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$',fontsize=18)
	ax = plt.subplot(gs[0,1])
	axes.append(ax)
	ax.set_title('term 2')
	im1 =ax.imshow(array(t1).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}C_{0,-q}\right\rangle}$',fontsize=18)
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(im1,cax=cax1)
	ax = plt.subplot(gs[1,0])
	axes.append(ax)
	ax.set_title('term 3')
	ax.imshow(array(t2).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle C_{0,q}h_{-q}\right\rangle}$',fontsize=18)
	ax = plt.subplot(gs[1,1])
	axes.append(ax)
	ax.set_title('term 4')
	im3 = ax.imshow(array(t3).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle C_{0,q}C_{0,-q}\right\rangle}$',fontsize=18)
	divider = make_axes_locatable(ax)
	cax3 = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(im3,cax=cax3)
	for a in range(len(axes)):
		ax = axes[a]
		if a != 0:
			ax.set_yticks([])
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])	
		if a in [0,2]:
			ax.set_ylabel(r'$\left|\mathbf{q_y}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		if a in [2,3]:
			ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm^{-1}})$',fontsize=14)
	plt.savefig(pickles+'fig-bilayer-couple-view-terms-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()

#---plots, compare 2D undulation spectra between bare and protein systems
if int(seq[3]) and barecompare:
	fig = plt.figure(figsize=(8,8))
	gs2 = gridspec.GridSpec(1,2,wspace=0.4)
	axes = []
	m,n = mset.griddims
	imdat0 = (abs(mean(mset.uqraw,axis=0))/double(m*n)**2)
	m,n = mset2.griddims
	imdat1 = (abs(mean(mset2.uqraw,axis=0))/double(m*n)**2)
	extrema = [min(imdat0.min(),imdat1.min()),max(imdat0.max(),imdat1.max())]
	ax = plt.subplot(gs2[0])
	axes.append(ax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$, protein',fontsize=18)
	im = ax.imshow(array(imdat0).T, extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,norm=mpl.colors.LogNorm(),vmin=extrema[0],vmax=extrema[1])
	ax = plt.subplot(gs2[1])
	axes.append(ax)
	im = ax.imshow(array(imdat1).T,extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,norm=mpl.colors.LogNorm(),vmin=extrema[0],vmax=extrema[1])
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$, bare',fontsize=18)
	for a in range(len(axes)):
		ax = axes[a]
		ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		if a in [0]:		
			ax.set_ylabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=14)
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm+1)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn+1)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])
	plt.savefig(pickles+'fig-bilayer-couple-view-comparehqhq-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()

#---plots, compare 2D undulation spectra between bare and protein systems scaled by q4
if int(seq[4]) and barecompare:
	fig = plt.figure(figsize=(8,8))
	gs2 = gridspec.GridSpec(1,2,wspace=0.4)
	axes = []
	m,n = mset.griddims
	imdat0 = t0*qmagst**4
	m,n = mset2.griddims
	imdat1 = t0_bare*qmagst**4
	extrema = [min(imdat0.min(),imdat1.min()),max(imdat0.max(),imdat1.max())]
	ax = plt.subplot(gs2[0])
	axes.append(ax)
	ax.set_title(\
		r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle {\left|\mathbf{q_y}\right|}^{4}}$, protein',
		fontsize=18)
	im = ax.imshow(array(imdat0).T, extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,vmin=extrema[0],vmax=extrema[1])
	ax = plt.subplot(gs2[1])
	axes.append(ax)
	im = ax.imshow(array(imdat1).T,extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,vmin=extrema[0],vmax=extrema[1])
	ax.set_title(\
		r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle {\left|\mathbf{q_y}\right|}^{4}}$, bare',
		fontsize=18)
	for a in range(len(axes)):
		ax = axes[a]
		ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm}^{-1})$',fontsize=14)
		if a in [0]:			
			ax.set_ylabel(r'$\left|\mathbf{q_y}\right|(\mathrm{nm}^{-1})$',fontsize=14)
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])
	plt.savefig(pickles+'fig-bilayer-couple-view-comparehqhqq4-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()

#---summary of 1D spectra
if int(seq[5]):
	fig = plt.figure(figsize=(12,6))
	gs3 = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
	#---undulation spectra
	ax = plt.subplot(gs3[0])
	#---from plotter_undulate
	spec1d = array([i for i in array(mset.undulate_spec1d) if i[0] != 0.])
	specfilter = array(filter(lambda x: x[0] >= mset.undulate_qmagfilter[0]
		and x[0] <= mset.undulate_qmagfilter[1],spec1d))
	[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
	area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
	kappa = 1/exp(az)/area
	print 'kappa = '+str(kappa)+' kBT'
	leftcom = [mean(log(specfilter[:,0])),mean(log(specfilter[:,1]))]
	az_enforced = leftcom[1]+4.*leftcom[0]
	kappa_enforced = 1./exp(az_enforced)/area
	print 'kappa_enforced = '+str(kappa_enforced)
	#---plot settings
	ax.plot([10**-3,10**3],[exp(az_enforced)*(i**-4) for i in [10**-3,10**3]],c='k',lw=2,alpha=0.5)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((spec1d[:,0].min()/4,spec1d[:,0].max()*2))
	ax.set_ylim((t3spec[:,1].min()/4,max(spec1d[:,1])*2))
	ax.grid(True)
	ax.set_ylabel(r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle \left(\mathrm{nm}^{2}\right)$',
		fontsize=fsaxlabel)
	#---comparison to control		
	if barecompare:
		ax.scatter(specbare[:,0],specbare[:,1],color=clrs[4],marker='o',s=20,
			label=r'$\left\langle h_{q}h_{-q}\right\rangle$, bare')
	if mdcompare:
		ax.scatter(specmd[:,0],specmd[:,1],color=clrs[5],marker='o',s=20,
			label=r'$\left\langle h_{q}h_{-q}\right\rangle$, MD')
	ax.scatter(t0spec[:,0],t0spec[:,1],color=clrs[0],marker='o',s=20,
		label=r'$\left\langle h_{q}h_{-q}\right\rangle$'+
		'\n'+r'$\boldsymbol{\kappa} = '+str('%3.1f'%(kappa_enforced))+'\:k_BT$')
	#---plotting extra terms
	ax.scatter(t1spec[:,0],t1spec[:,1],color=clrs[1],marker='+')
	ax.scatter(t2spec[:,0],t2spec[:,1],color=clrs[2],marker='x',
		label=r'$\left\langle C_{0,q}h_{-q}\right\rangle$')
	ax.scatter(t3spec[:,0],t3spec[:,1],color=clrs[3],marker='x',s=20,
		label=r'$\left\langle C_{0,q}C_{0,-q}\right\rangle$')
	ax.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
	h,l = ax.get_legend_handles_labels()
	ax.set_title('undulations')
	if barecompare:
		h = [h[1],h[0],h[2],h[3]]
		l = [l[1],l[0],l[2],l[3]]
	plt.legend(h,l,loc='lower left')
	plt.tick_params(labelsize=fsaxticks)
	#---plot energy terms on the right panel
	ax = plt.subplot(gs3[1:])
	reord0 = np.lexsort((termsumspec[:,1],termsumspec[:,0]))
	ax.scatter(termsumspec[reord0,0],termsumspec[reord0,1]+termsumspec_err[reord0,1],
		color='0.8',s=20,marker='o')
	ax.fill_between(termsumspec[reord0,0],termsumspec[reord0,1]-termsumspec_err[reord0,1],
		termsumspec[reord0,1]+termsumspec_err[reord0,1],color='0.8',alpha=1.)
	corrects = array([[termsumspec[i,0],scalefacbare*t0spec2[i,1]] 
		for i in range(len(termsumspec))])
	reord1 = np.lexsort((corrects[:,1],corrects[:,0]))
	if barecompare:
		reord2 = np.lexsort((energy_bare[:,1],energy_bare[:,0]))
		ax.plot(energy_bare[reord2,0],energy_bare[reord2,1],'-',lw=1.5,color=clrs[4],alpha=0.5)	
		ax.scatter(energy_bare[reord2,0],energy_bare[reord2,1],marker='o',s=20,color=clrs[4],
			label=r'$\left\langle h_{q}h_{-q}\right\rangle {\left|q\right|}^4$, bare',)
	if mdcompare:
		reord3 = np.lexsort((energy_md[:,1],energy_md[:,0]))
		ax.plot(energy_md[reord3,0],energy_md[reord3,1],'-',lw=1.5,color=clrs[5],alpha=0.5)	
		ax.scatter(energy_md[reord3,0],energy_md[reord3,1],marker='o',s=20,color=clrs[5],
			label=r'$\left\langle h_{q}h_{-q}\right\rangle {\left|q\right|}^4$, MD',)
	ax.plot(termsumspec[reord1,0],corrects[reord1,1],'o-',color=clrs[0],lw=2,markersize=5,alpha=0.5)
	ax.scatter(termsumspec[reord1,0],corrects[reord1,1],marker='o',s=20,color=clrs[0],
		label=r'$\left\langle h_{q}h_{-q}\right\rangle {\left|q\right|}^4$')
	ax.plot(termsumspec[reord0,0],termsumspec[reord0,1],'-',color='k',lw=2,alpha=0.5)
	ax.scatter(termsumspec[reord0,0],termsumspec[reord0,1],marker='o',color='k',s=20,label='energy')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.yaxis.tick_right()
	ax.set_ylabel(\
		r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
		fontsize=fsaxlabel)
	ax.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
	ax.yaxis.set_label_position("right")
	h,l = ax.get_legend_handles_labels()
	plt.legend(h[::-1],l[::-1],loc='lower right')
	ax.grid(True,which='both')
	ax.set_xlim((min([i for i in energy_bare[:,0] if i != 0.])/(3./2),max(energy_bare[:,0])*(3./2)))
	ax.set_ylim((10**-2,4*10**1))
	ax.set_title('energy')
	plt.tick_params(labelsize=fsaxticks)
	plt.savefig(pickles+'fig-bilayer-couple-view-spectra1d-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()
	
#---phase plots
if int(seq[6]):
	#---analyze 2D plots of phase angles between several simulations
	msetmd2 = unpickle(pickles+'pkl.structures.membrane-v550-stress.md.part0008.shifted.pkl')
	msetmd2.calculate_undulations()
	if 1:
		gs = gridspec.GridSpec(1,4)
		for m in range(4):
			ms = [mset,mset2,msetmd,msetmd2][m]
			ax = plt.subplot(gs[m])
			angles = mean(array([angle(ms.undulate_raw[i]) for i in range(len(ms.undulate_raw))]),axis=0)
			anglesfilt = numpy.fft.ifftshift(scipy.ndimage.filters.gaussian_filter(angles,10))
			anglesfilt = scipy.ndimage.filters.gaussian_filter(angles,2)
			ax.imshow(anglesfilt.T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm(),
				cmap=mpl.cm.jet,vmax=pi,vmin=10**-10)
			ax.imshow(-1*anglesfilt.T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm(),
				cmap=mpl.cm.jet_r,vmax=pi,vmin=10**-10)
		plt.show()

	#---plot phase angle on XY for different wavelengths, or different groups of wavelengths
	#---Nb this just shows that the DTMC model isn't dynamic, obviously
	if 0:
		ind = 0
		fig = plt.figure()
		ax = plt.subplot(111,aspect='equal')
		wid = 5
		clrs = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,1,2,3]] #---RBG
		for m in range(4):
			ms = [msetmd,msetmd2,mset,mset2][m]
			dat = array([[mean(angle(ms.undulate_raw[i])[cm:cm+wid,cn]),
				mean(angle(ms.undulate_raw[i])[cm,cn:cn+wid])]
				for i in range(len(ms.undulate_raw))])
			plt.scatter(dat[:,0],dat[:,1],c=clrs[m])
		plt.show()


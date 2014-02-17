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
	('meso','meso-anis','v2002-t2-anis-22-run1-rep-0-1500-2000',
		'/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,22,
		True,None,True,True,False,[4,2],True),
	('meso','meso-bare','v2002-t2-bare-22-run1-rep-0-1500-2000',
		'/home/rpb/worker/repo-membrane/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,22,
		False,None,True,True,False,[4,2],True),
	('md','cgmd-ENTHx4','v614',
		'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl',None,None,None,
		True,[0.05,5,5,0],True,True,False,None,True)]
analyses_ord = [2,0,1]
plot_reord = [1,0,2]

analyses = [analysis_descriptors[i] for i in analyses_ord]
if 'msets' not in globals(): msets = []
if 'mscs' not in globals(): mscs = []
if 'collect_c0s' not in globals(): collect_c0s = []

#---method
removeavg = False
fitlims = [4,2]
forcekappa = True
match_scales = [2,0]

#---method 
'''	   |---------load
       ||--------calc
       |||-------plot terms
       ||||------plot compare hqhq
       |||||-----plot compare hqhqq4
       ||||||----plot spectrum, 1D
       |||||||---plot phase angles '''
seq = '0100010'

#---settings
#---load colors in the same order as analyses_ord
clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,1,2,3,4,5,6,7,]]
cmap = mpl.cm.jet
showplots = True
sskip = 4 #---sets the xy ticks spacing  in units of lenscale on any 2D spectra

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
		mset.calculate_undulations(removeavg=removeavg,fitlims=fitlims,forcekappa=forcekappa)
		grid = mset.griddims
		m,n = grid[0]-1,grid[1]-1
		if c0s == []:
			#---if no curvature field, use zeros
			c0s = zeros((m,n))
		hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
		cqs = [fftwrap(mset.lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
		Lx,Ly = mean(mset.vecs,axis=0)[0:2]
		cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
		qmags = mset.lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
			for j in range(0,n)] for i in range(0,m)])
		hqsa = autocorr(hqs,direct=0,lenscale=mset.lenscale)
		hqsb = autocorr(hqs,direct=1,lenscale=mset.lenscale)
		center = [cm,cn]
		cqsa = autocorr(cqs,direct=0,lenscale=mset.lenscale)
		cqsb = autocorr(cqs,direct=1,lenscale=mset.lenscale)
		qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
		mt,nt = shape(hqsa[0])
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
		scalefac = mset.undulate_kappa*area
		#---compute terms which contribute to the elastic energy in 2D
		self.t2d[0] = mean((abs(array(hqsa)))**2,axis=0)
		self.t2d[1] = mean(abs(hqsa)*abs(cqsb),axis=0)
		self.t2d[2] = mean(abs(cqsa*hqsb),axis=0)
		self.t2d[3] = mean(abs(cqsa*cqsb),axis=0)
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
		self.tsum2d = self.t2d[0]*qmagst**4-self.t2d[1]*qmagst**2-self.t2d[2]*qmagst**2+self.t2d[3]
		self.tsum1d = array([[qmagst[i,j],scalefac*self.tsum2d[i,j]] for j in range(nt) for i in range(mt)])
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
		self.c0s = c0s

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
if int(seq[0]) or msets == []:
	lenscale = 1.0
	for a in analyses_ord:
		[simtype,shortname,testname,locate,start,end,nbase,hascurv,
			hypo,plot_ener_err,plotqe,removeavg,fitlims,forcekappa] = analysis_descriptors[a]
		if 'mset' in globals(): del mset
		mset = MembraneSet()
		if simtype == 'meso':
			c0sraw = array(mset.load_points_vtu(locate,extra_props='induced_cur',
				start=start,end=end,nbase=nbase,lenscale=lenscale))[:,0]
			mset.surfacer()
			c0s = mset.surfacer_general(c0sraw)
			collect_c0s.append(c0s)
			msets.append(mset)
		elif simtype == 'md':
			print 'loading from MD'
			if 'mset' in globals(): del mset
			mset = unpickle(pickles+locate)
			msets.append(mset)
			#---compute hypothetical curvature field for the MD simulation
			if hypo != None:
				vecs = mean(mset.vecs,axis=0)
				m,n = mset.griddims
				getgrid = array([[[i,j] for j in linspace(0,vecs[1],n)] for i in linspace(0,vecs[0],m)])
				params = [0,hypo[0],vecs[0]/2.,vecs[1]/2.,hypo[1],hypo[2],hypo[3]]
				c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1])
					for j in range(n)] for i in range(m)])
				collect_c0s.append([c0hypo for i in range(len(mset.surf))])
			else:
				collect_c0s.append([])

	'''
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
	#lenscale = max(mean(mset2.vecs,axis=0))/(max(mean(msetmd.vecs,axis=0))/msetmd.lenscale)
	#print 'reset lenscale = '+str(lenscale)
	#mset.lenscale = lenscale
	#if barecompare:
	#	mset2.lenscale = lenscale
	'''
		
#---calculate mode couplings
if int(seq[1]):
	#---match mesoscale length scales to an MD simulation
	if match_scales != None:
		lenscale = max(mean(msets[analyses_ord.index(match_scales[1])].vecs,axis=0))/\
			(max(mean(msets[analyses_ord.index(match_scales[0])].vecs,axis=0))/msets[analyses_ord.index(match_scales[0])].lenscale)
		for a in range(len(analyses_ord)):
			if analysis_descriptors[analyses_ord[a]][0] == 'meso':
				msets[a].lenscale = lenscale
	#---calculate coupled modes
	for m in range(len(msets)):
		if 'msc' in globals(): del msc
		msc = ModeCouple()
		msc.calculate_mode_coupling(msets[m],collect_c0s[m])
		mscs.append(msc)
		
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
	#---prepare figure
	fig = plt.figure(figsize=(12,6))
	gs3 = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
	axl = plt.subplot(gs3[0])
	axr = plt.subplot(gs3[1])
	axl_range = [[10**-10,10**10],[10**-10,10**10]]
	for a in [analyses_ord.index(i) for i in plot_reord]:
		[simtype,shortname,testname,locate,start,end,
			nbase,hascurv,hypo,plot_ener_err,plotqe,
			removeavg,fitlims,forcekappa] = analysis_descriptors[analyses_ord[a]]
		colord = clist[a]
		mset = msets[a]
		#---calculate fitted line
		spec1d = array([i for i in array(mset.undulate_spec1d) if i[0] != 0.])
		specfilter = array(filter(lambda x: x[0] >= mset.undulate_qmagfilter[0]
			and x[0] <= mset.undulate_qmagfilter[1],spec1d))
		[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
		kappa = 1/exp(az)/area
		leftcom = [mean(log(specfilter[:,0])),mean(log(specfilter[:,1]))]
		az_enforced = leftcom[1]+4.*leftcom[0]
		kappa_enforced = 1./exp(az_enforced)/area
		axl.plot([10**-3,10**3],[exp(az_enforced)*(i**-4) for i in [10**-3,10**3]],c='k',lw=2,alpha=0.5)
		#---plot the undulations
		axl.scatter(mscs[a].t1d[0][:,0],mscs[a].t1d[0][:,1],color=colord,marker='o',s=20,
			label=r'$\left\langle h_{q}h_{-q}\right\rangle$'+', '+shortname+
			'\n'+r'$\boldsymbol{\kappa} = '+str('%3.1f'%(kappa_enforced))+'\:k_BT$')
		if axl_range[0][0] > min(mscs[a].t1d[0][:,0]): axl_range[0][0] = min(mscs[a].t1d[0][:,0])
		if axl_range[0][1] < max(mscs[a].t1d[0][:,0]): axl_range[0][1] = max(mscs[a].t1d[0][:,0])
		if axl_range[1][0] > min(mscs[a].t1d[0][:,1]): axl_range[1][0] = min(mscs[a].t1d[0][:,1])
		if axl_range[1][1] < max(mscs[a].t1d[0][:,1]): axl_range[1][0] = max(mscs[a].t1d[0][:,1])
		#---plotting extra terms
		if hascurv:
			axl.scatter(mscs[a].t1d[1][:,0],mscs[a].t1d[1][:,1],color=colord,marker='+',alpha=0.65)
			label = r'$\left\langle C_{0,q}h_{-q}\right\rangle$'+','+shortname  if 0 else None
			axl.scatter(mscs[a].t1d[2][:,0],mscs[a].t1d[2][:,1],color=colord,marker='x',
				label=label,alpha=0.65)
			label = r'$\left\langle C_{0,q}C_{0,-q}\right\rangle$'+','+shortname if 0 else None
			axl.scatter(mscs[a].t1d[3][:,0],mscs[a].t1d[3][:,1],color=colord,marker='x',s=20,
				label=label,alpha=0.65)
		#---plot details
		axl.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
		h,l = axl.get_legend_handles_labels()
		axl.set_title('undulations')
		axl.legend(h[::-1],l[::-1],loc='lower left')
		axl.tick_params(labelsize=fsaxticks)	
		axl.set_xscale('log')
		axl.set_yscale('log')
		axl.grid(True)
		axl.set_ylabel(\
			r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle \left(\mathrm{nm}^{2}\right)$',
			fontsize=fsaxlabel)
		#---plot energies
		reord0 = np.lexsort((mscs[a].tsum1d[:,1],mscs[a].tsum1d[:,0]))
		if plot_ener_err:
			axr.fill_between(mscs[a].tsum1d[reord0,0],mscs[a].tsum1d[reord0,1]-mscs[a].tsum1de[reord0,1],
				mscs[a].tsum1d[reord0,1]+mscs[a].tsum1de[reord0,1],color=colord,alpha=0.2)
		#---plot the height-height autocorrelation-equivalent energy function on the right panel
		if plotqe:
			reord0b = np.lexsort((mscs[a].t1denergy[0][:,1],mscs[a].t1denergy[0][:,0]))
			axr.plot(mscs[a].t1denergy[0][reord0b,0],mscs[a].t1denergy[0][reord0b,1],'--',
				color=colord,lw=2,alpha=1.)
			axr.scatter(mscs[a].t1denergy[0][reord0b,0],mscs[a].t1denergy[0][reord0b,1],
				marker='o',color=colord,s=20)
		#---plot the data
		axr.plot(mscs[a].tsum1d[reord0,0],mscs[a].tsum1d[reord0,1],'-',color=colord,lw=2,alpha=0.5)
		axr.scatter(mscs[a].tsum1d[reord0,0],mscs[a].tsum1d[reord0,1],
			marker='o',color=colord,s=20,label=shortname)
		#---plot details
		axr.set_xscale('log')
		axr.set_yscale('log')
		axr.yaxis.tick_right()
		axr.set_ylabel(\
			r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
			fontsize=fsaxlabel)
		axr.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
		axr.yaxis.set_label_position("right")
		h,l = axr.get_legend_handles_labels()
		plt.legend(h[::-1],l[::-1],loc='lower right')
		axr.grid(True,which='both')
		axr.set_title('energy')
		plt.tick_params(labelsize=fsaxticks)
	#---set final plots limits for the undulations
	xlims = (
		1./2*min([min([min([i for i in mscs[analyses_ord.index(k)].t1d[t][:,0] if i != 0.]) 
			for t in ([0,1,2,3] 
			if analysis_descriptors[k][7] else [0])]) for k in analyses_ord]),
		2*max([max([max([i for i in mscs[analyses_ord.index(k)].t1d[t][:,0] if i != 0.])
			for t in ([0,1,2,3] 
			if analysis_descriptors[k][7] else [0])]) for k in analyses_ord])
		)
	axl.set_xlim(xlims)
	axr.set_xlim(xlims)
	axr.set_ylim((1./2*10**-1,2*10**1))
	axl.set_ylim((
		1./2*min([min([min([i for i in mscs[analyses_ord.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if analysis_descriptors[k][7] else [0])]) for k in analyses_ord]),
		2*max([max([max([i for i in mscs[analyses_ord.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if analysis_descriptors[k][7] else [0])]) for k in analyses_ord]),
		))
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


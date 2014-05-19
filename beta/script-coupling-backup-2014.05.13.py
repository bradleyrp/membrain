#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

import scipy.ndimage
from scipy.signal import hilbert

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
'''
analysis_descriptors = {
	'v2002-t3':
		{'simtype':'meso_precomp',
		'shortname':r'meso(iso)',
		'testname':'v2002-t3',
		'locate':'pkl.structures.meso.v2002-t3.pkl',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':True,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t4':
		{'simtype':'meso_precomp',
		'shortname':'meso(bare)',
		'testname':'v2002-t4',
		'locate':'pkl.structures.meso.v2002-t4.pkl',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':False,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t2':
		{'simtype':'meso_precomp',
		'shortname':r'meso(aniso)',
		'testname':'v2002-t2',
		'locate':'pkl.structures.meso.v2002-t2.pkl',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':True,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v2002-t1':
		{'simtype':'meso_precomp',
		'shortname':'meso(bare)',
		'testname':'v2002-t1',
		'locate':'pkl.structures.meso.v2002-t1.pkl',
		'start':1500,
		'end':2000,
		'nbase':22,
		'hascurv':False,
		'hypo':None,
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':[16,4],
		'forcekappa':True},
	'v614':
		{'simtype':'md',
		'shortname':r'$4\times$ENTH(MD)',
		'testname':'v614',
		'locate':'pkl.structures.membrane-v614.s9-lonestar.md.part0004.120000-220000-200.pkl',
		'start':None,
		'end':None,
		'nbase':None,
		'hascurv':True,
		'hypo':[0.005,1./2,1./2,5,5,0],
		'plot_ener_err':True,
		'plotqe':True,
		'removeavg':False,
		'fitlims':None,
		'forcekappa':True}}
'''

#---load the standard header defintions
execfile('header-meso.py')

analyses_names = ['v2004-C_0-0.005'][:]

#---analyses
dotest = 'v2004-C_0-0.005.v614'
if dotest == 'v2002.t3.t4.v614':
	analyses_names = ['v614','v2002-t3']
	#---Nb I think 'v2002-t3' is missing c0s in the mset or they are zeros everywhere
	#---reverse order of importance for the 1D spectrum plot
	plot_reord = ['v2002-t3','v614']
	match_scales = ['v614','v2002-t3']
	analysis_name = 'v2002.t3.v614'
	#---load colors in the same order as analyses_names
	#---previous v614 hypo [0.005,2./5,2./5,10,10,0]
	clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,2,1,3,4,5,6,7,]]
	routine = ['load','calc','checkplot','plot1d','plot2d','plotphase','hyposweep'][1:3]
	index_md = 0
	if 'hyposweep' in routine:
		hypo_list = [[i,1./2,1./2,10,10,0] for i in arange(0.1,1.0,0.05)]
		hypo_list = [[0.6,1./2,1./2,i,i,0] for i in range(20)+list(arange(30,100,10))]
		hypo_list = [[0.6,i,j,14,14,0] for i in arange(0,1,0.05) for j in arange(0,1,0.05)]
		hypo_list = [[0.6,i,j,14,14,0] for i in arange(0,1,0.05) for j in arange(0,1,0.05)]
		hypo_list = [[0.6,1./2,1./2,i*j,i,pi/4] for j in [1,2,4,8,10] for i in list(arange(2,20,2))]
		hypo_list = [[i,cx,cy,sig1,sig2,0] 
			for cx in [0.5] 
			for cy in [0.5] 
			for sig1 in range(1,10+1,2) 
			for sig2 in range(1,10+1,2)
			for i in arange(0.01,0.04,0.05)]
elif dotest == 'v2002.t2.t1':
	analyses_names = ['v2002-t2','v2002-t1']
	plot_reord = ['v2002-t1','v2002-t2']
	match_scales = None
	analysis_name = 'v2002.t2.t1'
	#---load colors in the same order as analyses_names
	clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,1,2,3,4,5,6,7,]]
elif dotest == 'v2004-C_0-0.005.v614':
	analyses_names = ['v614','v2004-C_0-0.005']
	#---Nb I think 'v2002-t3' is missing c0s in the mset or they are zeros everywhere
	#---reverse order of importance for the 1D spectrum plot
	plot_reord = ['v2004-C_0-0.005','v614']
	match_scales = ['v614','v2004-C_0-0.005']
	analysis_name = 'v2004-C_0-0.005.v614'
	#---load colors in the same order as analyses_names
	#---previous v614 hypo [0.005,2./5,2./5,10,10,0]
	clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,2,1,3,4,5,6,7,]]
	routine = ['load','calc','checkplot','plot1d','plot2d','plotphase','hyposweep'][1:3]
	index_md = 0
	if 'hyposweep' in routine:
		hypo_list = [[i,1./2,1./2,10,10,0] for i in arange(0.1,1.0,0.05)]
		hypo_list = [[0.6,1./2,1./2,i,i,0] for i in range(20)+list(arange(30,100,10))]
		hypo_list = [[0.6,i,j,14,14,0] for i in arange(0,1,0.05) for j in arange(0,1,0.05)]
		hypo_list = [[0.6,i,j,14,14,0] for i in arange(0,1,0.05) for j in arange(0,1,0.05)]
		hypo_list = [[0.6,1./2,1./2,i*j,i,pi/4] for j in [1,2,4,8,10] for i in list(arange(2,20,2))]
		hypo_list = [[i,cx,cy,sig1,sig2,0] 
			for cx in [0.5] 
			for cy in [0.5] 
			for sig1 in range(1,10+1,2) 
			for sig2 in range(1,10+1,2)
			for i in arange(0.01,0.04,0.05)]

if 'msets' not in globals(): msets = []; plt.ion()
if 'mscs' not in globals(): mscs = []
if 'collect_c0s' not in globals(): collect_c0s = []

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

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
		#---if no curvature field, use zeros
		if c0s == []: c0s = zeros((m,n))
		#---units of hqs are unspecified here but then set by autocorr
		hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
		#---units of incoming c0s must be in the units specified by mset.lenscale
		cqs = [fftwrap(mset.lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
		Lx,Ly = mean(mset.vecs,axis=0)[0:2]
		cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
		qmags = mset.lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
			for j in range(0,n)] for i in range(0,m)])
		hqsa = autocorr(hqs,direct=0,lenscale=mset.lenscale)
		hqsb = autocorr(hqs,direct=1,lenscale=mset.lenscale)
		center = [cm,cn]
		cqsa = autocorr(cqs,direct=0,lenscale=1.0)
		cqsb = autocorr(cqs,direct=1,lenscale=1.0)
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
		self.tsum2d = scalefac*(self.t2d[0]*qmagst**4-self.t2d[1]*qmagst**2-
			self.t2d[2]*qmagst**2+self.t2d[3])
		self.tsum1d = array([[qmagst[i,j],self.tsum2d[i,j]] for j in range(nt) for i in range(mt)])
		self.tsum2de = array([[sqrt((self.t2de[0][i,j]*qmagst[i,j]**4)**2+
			(qmagst[i,j]**2*self.t2de[1][i,j])**2+
			(qmagst[i,j]**2*self.t2de[2][i,j])**2+
			(self.t2de[3][i,j])**2) 
			for j in range(nt)] for i in range(mt)])
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

def spectrum_summary(hypotext=False):
	'''Plots the 1D spectrum from global variables.'''
	#---prepare figure
	fig = plt.figure(figsize=(18,8))
	gs3 = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
	axl = plt.subplot(gs3[0])
	axr = plt.subplot(gs3[1])
	axl_range = [[10**-10,10**10],[10**-10,10**10]]
	#---plot undulations
	for m in [analyses_names.index(aname) for aname	in plot_reord]:
		a = analyses_names[m]
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		colord = clist[m]
		mset = msets[m]
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
		axl.scatter(mscs[m].t1d[0][:,0],mscs[m].t1d[0][:,1],color=colord,marker='o',s=40,
			label=r'$\left\langle h_{q}h_{-q}\right\rangle$'+', '+shortname+
			'\n'+r'$\boldsymbol{\kappa} = '+str('%3.1f'%(kappa_enforced))+'\:k_BT$')
		if axl_range[0][0] > min(mscs[m].t1d[0][:,0]): axl_range[0][0] = min(mscs[m].t1d[0][:,0])
		if axl_range[0][1] < max(mscs[m].t1d[0][:,0]): axl_range[0][1] = max(mscs[m].t1d[0][:,0])
		if axl_range[1][0] > min(mscs[m].t1d[0][:,1]): axl_range[1][0] = min(mscs[m].t1d[0][:,1])
		if axl_range[1][1] < max(mscs[m].t1d[0][:,1]): axl_range[1][0] = max(mscs[m].t1d[0][:,1])
		#---plotting extra terms
		if hascurv:
			axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1],color=colord,marker='+',alpha=0.65)
			label = r'$\left\langle C_{0,q}h_{-q}\right\rangle$'+','+shortname  if 0 else ''
			axl.scatter(mscs[m].t1d[2][:,0],mscs[m].t1d[2][:,1],color=colord,marker='x',
				label=label,alpha=0.65)
			label = r'$\left\langle C_{0,q}C_{0,-q}\right\rangle$'+','+shortname if 0 else ''
			axl.scatter(mscs[m].t1d[3][:,0],mscs[m].t1d[3][:,1],color=colord,marker='.',s=20,
				label=label,alpha=0.65)
		#---plot details
		axl.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
		h,l = axl.get_legend_handles_labels()
		h = [h[i] for i in range(len(l)) if l[i] != '']
		l = [l[i] for i in range(len(l)) if l[i] != '']
		axl.set_title('undulations')
		axl.legend(h[::-1],l[::-1],loc='lower left')
		axl.tick_params(labelsize=fsaxticks)	
		axl.set_xscale('log')
		axl.set_yscale('log')
		axl.grid(True)
		axl.set_ylabel(\
			r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle \left(\mathrm{nm}^{2}\right)$',
			fontsize=fsaxlabel)
	#---plot energy errors
	for m in [analyses_names.index(aname) for aname	in plot_reord]:
		a = analyses_names[m]
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		colord = clist[m]
		mset = msets[m]
		reord0 = np.lexsort((mscs[m].tsum1d[:,1],mscs[m].tsum1d[:,0]))
		if plot_ener_err:
			if 0:
				axr.fill_between(mscs[m].tsum1d[reord0,0],
					mscs[m].tsum1d[reord0,1]-mscs[m].tsum1de[reord0,1],
					mscs[m].tsum1d[reord0,1]+mscs[m].tsum1de[reord0,1],color=colord,alpha=0.2)
			if 1:
				#---use hilbert to get the analytical signal to plot error bars more intelligently
				axr.fill_between(mscs[m].tsum1d[reord0[::1],0],
					real(hilbert(mscs[m].tsum1d[reord0[::-1],1]-mscs[m].tsum1de[reord0[::-1],1])),
					real(hilbert(mscs[m].tsum1d[reord0[::1],1]+mscs[m].tsum1de[reord0[::1],1])),
					color=colord,alpha=0.2)
			if 0:
				axr.fill_between(mscs[m].tsum1d[reord0[::1],0],
					-1*real(hilbert(-1*(mscs[m].tsum1d[reord0[::1],1]-mscs[m].tsum1de[reord0[::1],1])))[::1],
					real(hilbert(1*(mscs[m].tsum1d[reord0[::1],1]+mscs[m].tsum1de[reord0[::1],1]))),
					color=colord,alpha=0.2)
	for m in [analyses_names.index(aname) for aname	in plot_reord]:
		a = analyses_names[m]
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		colord = clist[m]
		mset = msets[m]
		reord0 = np.lexsort((mscs[m].tsum1d[:,1],mscs[m].tsum1d[:,0]))
		if plotqe:
			reord0b = np.lexsort((mscs[m].t1denergy[0][:,1],mscs[m].t1denergy[0][:,0]))
			if 0:
				axr.plot(mscs[m].t1denergy[0][reord0b,0],mscs[m].t1denergy[0][reord0b,1],'--',
					color=colord,lw=2,alpha=1.)
			axr.scatter(mscs[m].t1denergy[0][reord0b,0],mscs[m].t1denergy[0][reord0b,1],
				marker='o',color=colord,s=40)
		#---plot the data
		if 0:
			axr.plot(mscs[m].tsum1d[reord0,0],mscs[m].tsum1d[reord0,1],'-',color=colord,lw=2,alpha=0.5)
		axr.scatter(mscs[m].tsum1d[reord0,0],mscs[m].tsum1d[reord0,1],
			marker='o',color=colord,s=40,label=shortname)
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
		axr.set_title('energy',fontsize=fsaxlabel)
		plt.tick_params(labelsize=fsaxticks)
	#---set final plots limits for the undulations
	xlims = (1./2*min([min([min([i for i in mscs[analyses_names.index(k)].t1d[t][:,0] if i != 0.]) 
		for t in ([0,1,2,3] if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analyses_names]),
		2*max([max([max([i for i in mscs[analyses_names.index(k)].t1d[t][:,0] if i != 0.]) for t in ([0,1,2,3] 
		if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analyses_names]))
	axl.set_xlim(xlims)
	axr.set_xlim(xlims)
	axr.set_ylim((1./2*10**-1,2*10**1))
	axl.set_ylim((
		1./2*min([min([min([i for i in mscs[analyses_names.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analyses_names]),
		2*max([max([max([i for i in mscs[analyses_names.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analyses_names])))
	#---fix the limits for now
	if 1: axl.set_ylim((10**-11,10**0))
	hypo_suffix = '' if not hypo else '-hypo-'+str(hypo[0])+'-'+str(hypo[1])+'-'+str(hypo[2])+'-'+\
		str(hypo[3])+'-'+str(hypo[4])+'-'+str(hypo[5])
	if hypotext:
		axl.text(0.95,0.95,
			r'$C_0='+str(hypo[0])+'\:\mathrm{{nm}^{-1}}$\n'+\
			r'$(\sigma_{a},\sigma_{b})=('+str(hypo[3])+','+str(hypo[4])+')\:\mathrm{nm}$\n'+\
			r'$(\frac{x_0}{L_x},\frac{y_0}{L_y})=('+str(hypo[1])+','+str(hypo[2])+')$',
			horizontalalignment='right',
			verticalalignment='top',
			transform=axl.transAxes)
	plt.savefig(pickles+'fig-bilayer-couple-'+analysis_name+'-'+'spectrum1d'+hypo_suffix+'.png',
		dpi=200,bbox_inches='tight')
	plt.show()

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
if 'load' in routine or msets == []:
	lenscale = 1.0
	for a in analyses_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if 'mset' in globals(): del mset
		mset = MembraneSet()
		if simtype == 'meso':
			c0sraw = array(mset.load_points_vtu(locate,extra_props='induced_cur',start=start,end=end,nbase=nbase,lenscale=lenscale))[:,0]
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
				#---getgrid is xyz points in nm
				getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
					for i in linspace(0,vecs[0]/mset.lenscale,m)])
				#---convert everything to nm
				#---recall z0,c0,x0,y0,sx,sy,th = params
				#---params sets C0 in native units, x0,y0 in proportional units, and sx,sy in nm
				params = [0,hypo[0],vecs[0]*hypo[1]/mset.lenscale,vecs[1]*hypo[2]/mset.lenscale,
					hypo[3],hypo[4],hypo[5]]
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

#---calculate mode couplings
if 'calc' in routine:
	#---match mesoscale length scales to an MD simulation
	if match_scales != None:
		ref_ind = analyses_names.index(match_scales[0])
		move_ind = analyses_names.index(match_scales[1])
		lenscale = max(mean(msets[move_ind].vecs,axis=0))/\
			(max(mean(msets[ref_ind].vecs,axis=0))/msets[ref_ind].lenscale)
		for a in analyses_names:
			for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
			if (analysis_descriptors[a])['simtype'] == 'meso' or \
				(analysis_descriptors[a])['simtype'] == 'meso_precomp':
				msets[analyses_names.index(a)].lenscale = lenscale
		#---reset the hypothetical C0 field according to the new scaling (replaces the calculation above)
		for a in analyses_names:
			for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
			if analysis_descriptors[a]['simtype'] == 'md':
				anum = analyses_names.index(a)
				mset = msets[anum]
				vecs = mean(mset.vecs,axis=0)
				m,n = mset.griddims
				params = [0,0.5*hypo[0]*msets[1].lenscale,vecs[0]*hypo[1]/mset.lenscale,vecs[1]*hypo[2]/mset.lenscale,hypo[3]/msets[1].lenscale,hypo[4]/msets[1].lenscale,hypo[5]]
				c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
				collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
	#---calculate coupled modes
	for a in analyses_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analyses_names.index(a)
		if 'msc' in globals(): del msc
		msc = ModeCouple()
		msc.calculate_mode_coupling(msets[m],collect_c0s[m])
		mscs.append(msc)
	hypo = (analysis_descriptors['v614'])['hypo']
	spectrum_summary(hypotext=True)
	#execfile('script-coupling-dev2.py')
		
#---calculate mode couplings according to testable hypotheses
if 'hyposweep' in routine:
	#---match mesoscale length scales to an MD simulation
	if match_scales != None:
		ref_ind = analyses_names.index(match_scales[0])
		move_ind = analyses_names.index(match_scales[1])
		lenscale = max(mean(msets[move_ind].vecs,axis=0))/\
			(max(mean(msets[ref_ind].vecs,axis=0))/msets[ref_ind].lenscale)
		for a in analyses_names:
			if (analysis_descriptors[a])['simtype'] == 'meso' or \
				(analysis_descriptors[a])['simtype'] == 'meso_precomp':
				msets[analyses_names.index(a)].lenscale = lenscale
	#---calculate coupled modes
	for a in analyses_names:
		m = analyses_names.index(a)
		if 'msc' in globals(): del msc
		msc = ModeCouple()
		msc.calculate_mode_coupling(msets[m],collect_c0s[m])
		mscs.append(msc)
	#---loop over testable hypotheses
	mset = msets[index_md]
	for hypo in hypo_list:
		print 'status: testing hypothesis: '+str(hypo)
		vecs = mean(mset.vecs,axis=0)
		m,n = mset.griddims
		getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] for i in linspace(0,vecs[0]/mset.lenscale,m)])
		#---convert everything to nm
		#---recall z0,c0,x0,y0,sx,sy,th = params
		params = [0,hypo[0],vecs[0]*hypo[1]/mset.lenscale,vecs[1]*hypo[2]/mset.lenscale,hypo[3],hypo[4],hypo[5]]
		c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
		mscs[index_md].calculate_mode_coupling(msets[index_md],[c0hypo for i in range(len(mset.surf))])
		spectrum_summary(hypotext=True)		

#---plots, compare 2D undulation spectra between bare and protein systems alone, or scaled by q4
if 'plot2d' in routine:
	insets = True
	i2wid = 1
	#---plot these for both the <h_{q}h_{-q}> and <h_{q}h_{-q}>q4
	for d in ['helastic','hq2','helerr','hq2q4']:
		fig = plt.figure(figsize=(4*len(analyses_names),6))
		gs = gridspec.GridSpec(1,len(analyses_names),hspace=0.5,wspace=0.5)
		for m in [analyses_names.index(aname) for aname in plot_reord[::-1]]:
			a = analyses_names[m]
			mset = msets[m]
			#---note that you can also plot the errors by using mscs[m].tsum2de
			if d == 'helastic':
				data = mscs[m].tsum2d
				if (analysis_descriptors[a])['simtype'] == 'md':
					cm,cn = [int(i/2)-1 for i in shape(mset.undulate_hqhq2d)]
					data[cm,cn] = 1.0
				title = str((analysis_descriptors[a])['shortname'])+'\n'+\
					r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$'
				lims = [1*10**-1,1*10**1]
				cmap = mpl.cm.RdBu_r
			elif d == 'hq2': 
				data = mscs[m].t2d[0]
				title = str((analysis_descriptors[a])['shortname'])+'\n'+\
					r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$'
				lims = None
				cmap = mpl.cm.jet
			elif d == 'helerr':
				data = mscs[m].tsum2de
				if (analysis_descriptors[a])['simtype'] == 'md':
					cm,cn = [int(i/2)-1 for i in shape(mset.undulate_hqhq2d)]
					data[cm,cn] = mean(data)
				title = str((analysis_descriptors[a])['shortname'])+'\n'+\
					r'$\mathbf{\delta}(\mathbf{\left\langle h_{\mathbf{q}}'+\
					r'h_{\mathbf{\mathbf{-q}}}\right\rangle q}^{4})$'
				lims = None
				cmap = mpl.cm.jet
			elif d == 'hq2q4':
				data = mscs[m].t1denergy[0]
				data = mscs[m].scalefac*(mscs[m].t2d[0]*mscs[m].qmagst**4)
				title = str((analysis_descriptors[a])['shortname'])+'\n'+\
					r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle {\left|\mathbf{q_y}\right|}^{4}}$'
				lims = None
				lims = [1*10**-1,1*10**1]
				cmap = mpl.cm.RdBu_r
			ax = plt.subplot(gs[(plot_reord[::-1]).index(a)])
			plotter2d(ax,mset,dat=data,cmap=cmap,lims=lims)
			#---Nb you can replace the following block with a single call to recursive function above
			if insets:
				cm,cn = [int(i/2)-1 for i in shape(mset.undulate_hqhq2d)]
				axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="30%",height="30%",loc=1)
				plotter2d(axins,mset,
					dat=data[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1],
					tickshow=False,cmap=cmap,lims=[array([i for i in flatten(data) 
						if i != 0.]).min(),data.max()])
			#ax.set_title(title,fontsize=fsaxlabel)
		plt.savefig(pickles+'fig-bilayer-couple-'+analysis_name+'-'+d+'.png',
			dpi=500,bbox_inches='tight')
		plt.close(fig)
	
#---phase plots
if 'plotphase' in routine:
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

#-------------------!!!!!!!!!!!!@@@@@@@@@@@@#########################
if 'checkplot' in routine:		
	#---calculations
	mset = msets[0]
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	#---getgrid is xyz points in nm
	if 0:
		getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] for i in linspace(0,vecs[0]/mset.lenscale,m)])
		if 'tmp' not in globals(): tmp = array(msets[0].surf)
		msets[0].surf = []
		for i in range(len(tmp)):
			msets[index_md].surf.append(1*tmp[i])
		maxpos = unravel_index((mean(msets[0].surf,axis=0)/msets[0].lenscale).argmax(),msets[0].griddims)
		maxposgrid = [maxpos[i]/vecs[i]*msets[0].griddims[i] for i in range(2)]
		#---hypothesis has C0 in native units
		#hypo = [0.005,0.5,0.5,5/msets[1].lenscale,5/msets[1].lenscale,0]
		#params = [0,hypo[0],vecs[0]*hypo[1]/mset.lenscale,vecs[1]*hypo[2]/mset.lenscale,hypo[3],hypo[4],hypo[5]]
		#c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
		#mscs[index_md].calculate_mode_coupling(msets[index_md],[c0hypo for i in range(len(mset.surf))])
	#---figure
	fig = plt.figure()
	gs = gridspec.GridSpec(3,4,wspace=0.0,hspace=0.0)
	#---plot curvature field
	extrem = max([max([j.max() for j in mscs[i].c0s]) for i in range(2)])
	ax = plt.subplot(gs[0,0])
	ax.imshow(mean(mscs[0].c0s,axis=0),vmax=extrem,vmin=0.,cmap=mpl.cm.binary)
	ax.set_title('CGMD')
	ax = plt.subplot(gs[0,1])
	ax.set_title('MESO')
	im = ax.imshow(mean(mscs[1].c0s,axis=0),vmax=extrem,vmin=0.,cmap=mpl.cm.binary)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels())
	plt.setp(axins.get_xticklabels())
	axins.set_ylabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',
		rotation=270)
	#---plot average structure
	vmax = max([mean(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	vmin = min([mean(msets[i].surf,axis=0).min()/msets[i].lenscale for i in range(2)])
	extrem = max(abs(vmax),abs(vmin))
	vmax,vmin = extrem,-extrem
	ax = plt.subplot(gs[1,0])
	ax.imshow(mean(msets[0].surf,axis=0)/msets[0].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r)
	ax = plt.subplot(gs[1,1])
	im = ax.imshow(mean(msets[1].surf,axis=0)/msets[1].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels())
	plt.setp(axins.get_xticklabels())
	axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle (\mathrm{nm})$',rotation=270)
	#---plot standard deviations
	extrem = max([std(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	ax = plt.subplot(gs[2,0])
	ax.imshow(std(msets[0].surf,axis=0)/msets[0].lenscale,vmin=0.,vmax=extrem,cmap=mpl.cm.jet)
	ax = plt.subplot(gs[2,1])
	im = ax.imshow(std(msets[1].surf,axis=0)/msets[1].lenscale,vmin=0.,vmax=extrem,cmap=mpl.cm.jet)
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	plt.setp(axins.get_yticklabels())
	plt.setp(axins.get_xticklabels())
	axins.set_ylabel(r'$\left\langle \left(z-\overline{z}\right)^{2} \right\rangle (\mathrm{{nm}^2})$',
		rotation=270)
	#---spectrum plot
	axl = plt.subplot(gs[0,3])
	m = 1
	axl = plt.subplot(133)
	axl.set_title(r'$\left\langle C_{0,q} h_{-q} \right\rangle $')
	axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1],marker='.',color='r',s=40,label='MESO')
	#axl.scatter(mscs[m].t1d[0][:,0],mscs[m].t1d[0][:,1]/msets[m].lenscale,marker='.',color='r',s=40,label='MESO')
	m = 0
	axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1],marker='.',color='b',s=40,label='CGMD')
	#axl.scatter(mscs[m].t1d[0][:,0],mscs[m].t1d[0][:,1]/msets[m].lenscale,marker='.',color='b',s=40,label='CGMD')
	axl.set_ylim((10**-11,10**-1))
	axl.set_xlim((0.06,1.5))
	axl.set_yscale('log')
	axl.set_xscale('log')
	axl.grid(True)
	axl.yaxis.set_ticks_position("right")
	axl.legend(loc='lower left')
	gs.tight_layout(fig,h_pad=0.,w_pad=1.0)
	plt.savefig(pickles+'fig-bilayer-couple-c0qhq-compare.png',bbox_inches='tight')
	plt.show()


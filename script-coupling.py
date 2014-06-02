#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')
execfile('header-meso.py')

import scipy.ndimage

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

#---if script is run in standalone mode
if 'batch_override' not in globals() or not batch_override:

	cgmd_avail = [
		'v614-120000-220000-200',
		'v616-210000-310000-200',
		]
		
	meso_avail = [
		'v2004',
		'v2005',
		'v2006',
		'v2008',
		]

	routine = [
		'calc',
		'masterplot',
		'checkplot',
		'plot2d',
		'plotphase'
		][:2]

	#---select a single cgmd experiment and a panel of meso experiments, from which c0ask will be used
	batch_cgmd = cgmd_avail[-1]
	batch_meso = meso_avail[3]
	c0ask = 0.0541
	
	#---set curvature extent from fixed parameter in the toc, assuming isotropic
	r_2 = (meso_expt_toc[batch_meso])['R_2']
	
	showplots = True
	
	#---bookkeeping
	p_interest = (meso_expt_toc[key])['parameter_name']
	analysis_names = [batch_cgmd,batch_meso+'-'+p_interest+'-'+str(c0ask)]
	plot_reord = analysis_names
	match_scales = [batch_cgmd,batch_meso+'-'+p_interest+'-'+str(c0ask)]
	bigname = '-'.join(analysis_names)

#---allocate if empty
if 'msets' not in globals(): msets = []; 
if showplots and 'msets' not in globals(): plt.ion()
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
		mset.calculate_undulations(removeavg=removeavg,fitlims=fitlims,
			forcekappa=forcekappa,qmagfilter=qmagfilter)
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

def spectrum_summary(fig=None,gs=None,lowlim=10**-14,titletext='undulations'):
	'''Plots the 1D spectrum from global variables.'''
	#---prepare figure
	savefig = False if fig != None else True
	if fig == None: fig = plt.figure(figsize=(18,8))
	if gs == None: gs3 = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
	else: gs3 = gs
	axl = plt.subplot(gs3[0])
	axr = plt.subplot(gs3[1])
	axl_range = [[10**-10,10**10],[10**-10,10**10]]
	clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [1,2,0,3,4,5,6,7,]]
	#---plot undulations
	extra_legend = []
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		a = analysis_names[m]
		shortname = (analysis_descriptors[a])['detail_name']
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
		axl.plot([10**-3,10**3],[exp(az)*(i**bz) for i in [10**-3,10**3]],'--',c=clist[m],
			lw=2,zorder=0)
		axl.plot([10**-3,10**3],[exp(az_enforced)*(i**-4) for i in [10**-3,10**3]],c='k',
			lw=2,zorder=0)
		#---plot the undulations
		if 0: xdat,ydat = mscs[m].t1d[0][:,0],mscs[m].t1d[0][:,1]
		#---new method
		xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
		ydat = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
		plotobj = axl.scatter(xdat,ydat,color=colord,marker='o',s=30,
			label=r'$\left\langle h_{q}h_{-q}\right\rangle$'+', '+shortname+
			'\n'+r'$\boldsymbol{\kappa} = '+str('%3.1f'%(kappa_enforced))+'\:k_BT$',
			zorder=1)
		axl.scatter([xdat[i] for i in range(len(xdat)) if xdat[i]<qmagfilter[1]],
			[ydat[i] for i in range(len(ydat)) if xdat[i]<qmagfilter[1]],color=colord,marker='o',s=30,
			zorder=2,lw=0.5,edgecolor='k')
		if m == 0: extra_legend.append(plotobj)
		if axl_range[0][0] > min(mscs[m].t1d[0][:,0]): axl_range[0][0] = min(mscs[m].t1d[0][:,0])
		if axl_range[1][0] > min(mscs[m].t1d[0][:,1]): axl_range[1][0] = min(mscs[m].t1d[0][:,1])
		if axl_range[0][1] < max(mscs[m].t1d[0][:,0]): axl_range[0][1] = max(mscs[m].t1d[0][:,0])
		if axl_range[1][1] < max(mscs[m].t1d[0][:,1]): axl_range[1][0] = max(mscs[m].t1d[0][:,1])
		#---plotting extra terms
		if hascurv:

			if 0: xdat,ydat = mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1]
			xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[1])

			plotobj = axl.scatter(xdat,ydat,facecolor=colord,marker='o',s=40,
				edgecolor='k',lw=0.5,alpha=0.65,zorder=2)
			if m == 0: extra_legend.append(plotobj)

			if 0: xdat,ydat = mscs[m].t1d[2][:,0],mscs[m].t1d[2][:,1]
			xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[2])

			axl.scatter(xdat,ydat,facecolor=colord,marker='o',s=40,
				edgecolor='k',lw=0.5,label='',alpha=0.65,zorder=3)

			if 0: xdat,ydat = mscs[m].t1d[3][:,0],mscs[m].t1d[3][:,1]
			xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[3])

			plotobj = axl.scatter(xdat,ydat,color=colord,marker=('x' if m==0 else '+'),s=20,
				label='',alpha=1,zorder=4,lw=1.5)
			if m == 0: extra_legend.append(plotobj)
		#---extra legend
		legend2 = axl.legend(extra_legend,
			[r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle$',
			r'$\left\langle C_{0,\mathbf{q}} h_{-\mathbf{q}} \right\rangle $',
			r'$\left\langle C_{0,\mathbf{q}} C_{0,-\mathbf{q}} \right\rangle $'],
			loc='upper right')
		for legobj in legend2.legendHandles:
		    legobj.set_color('k')
		    legobj.set_alpha(0.65)
		#---plot details
		axl.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
		h,l = axl.get_legend_handles_labels() 
		h = [h[i] for i in range(len(l)) if l[i] != '']
		l = [l[i] for i in range(len(l)) if l[i] != '']
		axl.set_title(titletext,fontsize=fsaxtitle)
		axl.legend(h[::-1],l[::-1],loc='lower left')
		plt.gca().add_artist(legend2)
		axl.tick_params(labelsize=fsaxticks)	
		axl.set_xscale('log')
		axl.set_yscale('log')
		axl.grid(True)
		axl.set_ylabel(\
			r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle \left(\mathrm{nm}^{2}\right)$',
			fontsize=fsaxlabel)
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		a = analysis_names[m]
		shortname = (analysis_descriptors[a])['detail_name']
		colord = clist[m]
		mset = msets[m]
		reord0 = np.lexsort((mscs[m].tsum1d[:,1],mscs[m].tsum1d[:,0]))
		if plotqe and 0:
			reord0b = np.lexsort((mscs[m].t1denergy[0][:,1],mscs[m].t1denergy[0][:,0]))
			axr.scatter(mscs[m].t1denergy[0][reord0b,0],mscs[m].t1denergy[0][reord0b,1],
				marker='o',color=colord,s=20,zorder=6)
		#---plot the data
		if 0: xdat,ydat = mscs[m].tsum1d[reord0,0],mscs[m].tsum1d[reord0,1]
		xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
		ydat = collapse_spectrum(mscs[m].qmagst,mscs[m].tsum2d)

		if 0: axr.scatter(xdat,ydat,
			marker='o',color=colord,s=20,label=shortname,zorder=7)
		axr.plot(xdat,ydat,'o-',color=colord,label=shortname)
		#---plot a single line
		if 0:
			inds = unique(mscs[m].tsum1d[:,0],return_inverse=True)[1]
			dat = mscs[m].tsum1d[:,0]
			xdat = [mean(dat[where(inds==i)[0]]) for i in range(max(inds))]
			dat = mscs[m].tsum1d[:,1]
			ydat = [mean(dat[where(inds==i)[0]]) for i in range(max(inds))]
			axr.plot(xdat,ydat,'o-',color=colord,label='')
		#---plot details
		axr.set_xscale('log')
		axr.set_yscale('log')
		axr.yaxis.tick_right()
		axr.set_ylabel(\
			r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
			fontsize=fsaxlabel,rotation=270)
		axr.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm}^{-1})$',fontsize=fsaxlabel)
		axr.yaxis.set_label_position("right")
		h,l = axr.get_legend_handles_labels()
		h = [h[i] for i in range(len(l)) if l[i] != '']
		l = [l[i] for i in range(len(l)) if l[i] != '']
		plt.legend(h[::-1],l[::-1],loc='lower right')
		axr.grid(True,which='both')
		axr.set_title('energy',fontsize=fsaxtitle)
		plt.tick_params(labelsize=fsaxticks)
	axl.grid(b=True,which='minor',color='k',linestyle=':',alpha=0.5)
	axl.axvline(x=qmagfilter[1],ymin=0,ymax=1.,lw=2,c='r',alpha=0.5,zorder=0)
	axr.axvline(x=qmagfilter[1],ymin=0,ymax=1.,lw=2,c='r',alpha=0.5,zorder=0)
	#---set final plots limits for the undulations
	xlims = (1./2*min([min([min([i for i in mscs[analysis_names.index(k)].t1d[t][:,0] if i != 0.]) 
		for t in ([0,1,2,3] if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analysis_names]),
		2*max([max([max([i for i in mscs[analysis_names.index(k)].t1d[t][:,0] if i != 0.]) 
			for t in ([0,1,2,3] 
		if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analysis_names]))
	axl.set_xlim(xlims)
	axr.set_xlim(xlims)
	axr.set_ylim((1./2*10**-1,2*10**1))
	axr.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
	axl.set_ylim((
		1./2*min([min([min([i for i in mscs[analysis_names.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analysis_names]),
		2*max([max([max([i for i in mscs[analysis_names.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if (analysis_descriptors[k])['hascurv'] else [0])]) for k in analysis_names])))
	#---fix the plot limits
	if 1: axl.set_ylim((lowlim,10**0))
	if savefig: plt.savefig(pickles+'fig-bilayer-couple-spectrum-'+bigname+'.png',
		dpi=300,bbox_inches='tight')
	if showplots: plt.show()
	
def collapse_spectrum(qs,hs):
	'''For rectangular systems, treat lateral dimensions as equivalent and take their average.'''
	cen = array([i/2 for i in shape(qs)])
	discdists = [[sum(array([i-cen[0],j-cen[1]])**2) 
		for j in range(shape(qs)[1])] for i in range(shape(qs)[0])]
	uniqs = unique(flatten(discdists),return_inverse=True)[1]
	collapsed = [mean(ravel(hs)[where(uniqs==i)]) for i in range(uniqs.max()) 
		if mean(ravel(qs)[where(uniqs==i)]) != 0]
	return collapsed
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
if ('load' in routine or ('msets' in globals() and msets == [])) and 'analysis_names' in globals():
	lenscale = 1.0
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if 'mset' in globals(): del mset
		mset = MembraneSet()
		if simtype == 'meso':
			mset = unpickle(pickles+'pkl.structures.meso.'+(analysis_descriptors[a])['testname']+'.pkl')
			c0s = mset.getdata('c0map').data
			collect_c0s.append(c0s)
			msets.append(mset)
		elif simtype == 'md':
			status('status: loading from MD')
			if 'mset' in globals(): del mset
			mset = unpickle(pickles+locate)
			msets.append(mset)
			#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
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

#---calculate mode couplings
if 'calc' in routine:
	#---match mesoscale length scales to an MD simulation
	if match_scales != None:
		ref_ind = analysis_names.index(match_scales[0])
		move_ind = analysis_names.index(match_scales[1])
		#---matching the average box vectors here by average and not maximum
		lenscale = mean(mean(msets[move_ind].vecs,axis=0)[:2])/\
			(mean(mean(msets[ref_ind].vecs,axis=0)[:2])/msets[ref_ind].lenscale)
		for a in analysis_names:
			for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
			if (analysis_descriptors[a])['simtype'] == 'meso' or \
				(analysis_descriptors[a])['simtype'] == 'meso_precomp':
				msets[analysis_names.index(a)].lenscale = lenscale
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
				#---needs checked, "key step used to have a 0.5*hypo[0] here possible due to convention"
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
	#---calculate coupled modes
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		if 'msc' in globals(): del msc
		msc = ModeCouple()
		msc.calculate_mode_coupling(msets[m],collect_c0s[m])
		mscs.append(msc)
	hypo = (analysis_descriptors[batch_cgmd])['hypo']
	hypo[0] = c0ask
	if 'masterplot' not in routine and 'simple_summary' in routine: spectrum_summary()
	if 'batch_override' in globals() and batch_override:
		if analysis_names[0]+'-'+p_interest+'-'+str(c0ask) not in master_spectrum_dict.keys():
			master_spectrum_dict[analysis_names[0]+'-'+p_interest+'-'+str(c0ask)] = {\
				'c0ask':c0ask,
				'cgmd_qs':mscs[0].qmagst,
				'cgmd_t2d':mscs[0].t2d,
				'lenscale':msets[0].lenscale,
				}
		if analysis_names[1] not in master_spectrum_dict.keys():
			master_spectrum_dict[analysis_names[1]] = {\
				'c0ask':c0ask,
				'meso_qs':mscs[1].qmagst,
				'meso_t2d':mscs[1].t2d,
				'lenscale':msets[1].lenscale,
				}
		
#---plots, compare 2D undulation spectra between bare and protein systems alone, or scaled by q4
if 'plot2d' in routine:
	insets = True
	i2wid = 1
	#---plot these for both the <h_{q}h_{-q}> and <h_{q}h_{-q}>q4
	for d in ['helastic','hq2','helerr','hq2q4']:
		fig = plt.figure(figsize=(4*len(analysis_names),6))
		gs = gridspec.GridSpec(1,len(analysis_names),hspace=0.5,wspace=0.5)
		for m in [analysis_names.index(aname) for aname in plot_reord[::-1]]:
			a = analysis_names[m]
			mset = msets[m]
			#---note that you can also plot the errors by using mscs[m].tsum2de
			if d == 'helastic':
				data = mscs[m].tsum2d
				if (analysis_descriptors[a])['simtype'] == 'md':
					cm,cn = [int(i/2)-1 for i in shape(mset.undulate_hqhq2d)]
					data[cm,cn] = 1.0
				title = str((analysis_descriptors[a])['detail_name'])+'\n'+\
					r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$'
				lims = [1*10**-1,1*10**1]
				cmap = mpl.cm.RdBu_r
			elif d == 'hq2': 
				data = mscs[m].t2d[0]
				title = str((analysis_descriptors[a])['detail_name'])+'\n'+\
					r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$'
				lims = None
				cmap = mpl.cm.jet
			elif d == 'helerr':
				data = mscs[m].tsum2de
				if (analysis_descriptors[a])['simtype'] == 'md':
					cm,cn = [int(i/2)-1 for i in shape(mset.undulate_hqhq2d)]
					data[cm,cn] = mean(data)
				title = str((analysis_descriptors[a])['detail_name'])+'\n'+\
					r'$\mathbf{\delta}(\mathbf{\left\langle h_{\mathbf{q}}'+\
					r'h_{\mathbf{\mathbf{-q}}}\right\rangle q}^{4})$'
				lims = None
				cmap = mpl.cm.jet
			elif d == 'hq2q4':
				data = mscs[m].t1denergy[0]
				data = mscs[m].scalefac*(mscs[m].t2d[0]*mscs[m].qmagst**4)
				title = str((analysis_descriptors[a])['detail_name'])+'\n'+\
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
		plt.savefig(pickles+'fig-bilayer-couple-'+bigname+'-'+d+'.png',
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
		if showplots: plt.show()
	#---plot phase angle on XY for different wavelengths, or different groups of wavelengths
	#---Nb this just shows that the DTMC model isn't dynamic, obviously
	if 0:
		ind = 0
		fig = plt.figure()
		ax = plt.subplot(111,aspect='equal')
		wid = 5
		clrs = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [0,1,2,3]]
		for m in range(4):
			ms = [msetmd,msetmd2,mset,mset2][m]
			dat = array([[mean(angle(ms.undulate_raw[i])[cm:cm+wid,cn]),
				mean(angle(ms.undulate_raw[i])[cm,cn:cn+wid])]
				for i in range(len(ms.undulate_raw))])
			plt.scatter(dat[:,0],dat[:,1],c=clrs[m])
		if showplots: plt.show()

#---comparison of curvature between MESO and CGMD methods
#---plots height-curvature correlation alongsize structure and variations
if 'checkplot' in routine:		
	#---fixed limits for making smooth gifs
	extremz_checkplot_global = -2.5,2.5
	extrems_checkplot_global = 2.0
	extremc0_checkplot_global = 0.02
	#---calculations
	mset = msets[0]
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	#---figure
	fig = plt.figure(figsize=(12,12))
	gs = gridspec.GridSpec(5,2,hspace=0.3)
	gs.update(left=0.0,right=0.45)
	gs2 = gridspec.GridSpec(5,3)
	gs2.update(left=0.5,right=1.0)
	axeslist = []
	#---plot curvature field
	extrem = max([max([j.max() for j in mscs[i].c0s]) for i in range(2)])
	if extremc0_checkplot_global != None: extrem = extremc0_checkplot_global
	ax = plt.subplot(gs[0,0])
	ax.imshow(mean(mscs[0].c0s,axis=0).T,vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
		interpolation='nearest',origin='lower')
	axeslist.append(ax)
	ax.set_title('CGMD')
	ax = plt.subplot(gs[0,1])
	axeslist.append(ax)
	ax.set_title('MESO')
	im = ax.imshow(mean(mscs[1].c0s,axis=0).T,vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',
		rotation=270)
	#---plot average structure
	vmax = max([mean(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	vmin = min([mean(msets[i].surf,axis=0).min()/msets[i].lenscale for i in range(2)])
	extrem = max(abs(vmax),abs(vmin))
	vmax,vmin = extrem,-extrem
	if extremz_checkplot_global != None: vmin,vmax = extremz_checkplot_global
	ax = plt.subplot(gs[1,0])
	axeslist.append(ax)
	im = ax.imshow(mean(msets[0].surf,axis=0).T/msets[0].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	ax = plt.subplot(gs[1,1])
	axeslist.append(ax)
	im = ax.imshow(mean(msets[1].surf,axis=0).T/msets[1].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle (\mathrm{nm})$',rotation=270)
	#axins.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=6))
	#---plot standard deviations
	extrem = max([std(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	if extrems_checkplot_global != None: extrem = extrems_checkplot_global
	ax = plt.subplot(gs[2,0])
	axeslist.append(ax)
	ax.imshow((std(msets[0].surf,axis=0)/msets[0].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	ax = plt.subplot(gs[2,1])
	axeslist.append(ax)
	im = ax.imshow((std(msets[1].surf,axis=0)/msets[1].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle \left(z-\overline{z}\right)^{2} \right\rangle (\mathrm{{nm}^2})$',
		rotation=270)
	#---compare 2D energy spectra
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		mset = msets[m]
		axr2 = plt.subplot(gs[3,m])
		cm,cn = [int(i/2) for i in shape(mscs[m].tsum2d)]
		wid = 3
		dat = mscs[m].tsum2d[cm-wid:cm+wid+1,cn-wid:cn+wid+1]
		dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
		im = plotter2d(axr2,mset,dat=dat,
			cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
			ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=10,label_style='q',lims=[0.1,10],
			tickskip=int(round(mset.griddims[0]/6,-1)))
	axins = inset_axes(axr2,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=axr2.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(
		r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
		fontsize=fsaxlabel,rotation=270)
	#---compare height-curvature correlation in 2D
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		mset = msets[m]
		axr3 = plt.subplot(gs[4,m])
		cm,cn = [int(i/2) for i in shape(mscs[m].t2d[1])]
		wid = 3
		dat = mscs[m].t2d[1][cm-wid:cm+wid+1,cn-wid:cn+wid+1]
		dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
		im = plotter2d(axr3,mset,dat=dat,
			cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
			ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=10,label_style='q',lims=[10**-5,10**-2],
			tickskip=int(round(mset.griddims[0]/6,-1)))
	axins = inset_axes(axr3,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=axr3.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(
		r'$\left\langle C_{0,q} h_{-q} \right\rangle $',
		fontsize=fsaxlabel,rotation=270)
	axeslist.append(axins)
	for ax in axeslist:
		plt.setp(ax.get_yticklabels(),fontsize=10)
		plt.setp(ax.get_xticklabels(),fontsize=10)	
	#---spectrum plot
	axl = plt.subplot(gs2[1:4,:])
	axl.set_ylabel(r'$\left\langle C_{0,q} h_{-q} \right\rangle $',fontsize=fsaxlabel)
	axl.set_title(r'$\mathrm{C_{0,hypo}}='+('{0:.3f}'.format(c0ask))+'a_0^{-1}'+\
		'='+('{0:.3f}'.format(c0ask*msets[1].lenscale))+'\:\mathrm{({nm}^{-1})}$')
	#---hard-coded requirement that the CGMD simulation is first on the list
	index_md = 0
	clist = [(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)[i] for i in [1,2,0,3,4,5,6,7,]]
	for m in range(1,len(analysis_names)):
		axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1],marker='o',color=clist[m],s=40,label='MESO')
	m = index_md
	axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1],marker='o',color=clist[m],s=40,label='CGMD')
	axl.set_ylim((10**-11,10**-1))
	axl.set_xlim((0.06,1.5))
	axl.set_yscale('log')
	axl.set_xscale('log')
	axl.grid(True)
	axl.yaxis.set_ticks_position("right")
	axl.legend(loc='lower left')
	#---save
	plt.savefig(pickles+'fig-bilayer-couple-compare-'+bigname+'.png',bbox_inches='tight',dpi=300)
	if showplots: plt.show()
	plt.close(fig)

#---comparison of curvature between MESO and CGMD methods
#---plots height-curvature correlation alongsize structure and variations
if 'masterplot' in routine:		
	#---fixed limits for making smooth gifs
	extremz_checkplot_global = -2.5,2.5
	extrems_checkplot_global = 2.0
	extremc0_checkplot_global = 0.02
	#---calculations
	mset = msets[0]
	vecs = mean(mset.vecs,axis=0)
	m,n = mset.griddims
	#---figure
	fig = plt.figure(figsize=(18,12))
	gs = gridspec.GridSpec(5,2,hspace=0.15)
	gs.update(left=0.0,right=0.3)
	gs2 = gridspec.GridSpec(1,2)
	gs2.update(left=0.4,right=1.0)
	axeslist = []
	#---plot curvature field
	extrem = max([max([j.max() for j in mscs[i].c0s]) for i in range(2)])
	if extremc0_checkplot_global != None: extrem = extremc0_checkplot_global
	ax = plt.subplot(gs[0,0])
	ax.imshow(mean(mscs[0].c0s,axis=0).T,vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
		interpolation='nearest',origin='lower')
	axeslist.append(ax)
	ax.set_title('CGMD')
	ax = plt.subplot(gs[0,1])
	axeslist.append(ax)
	ax.set_title('MESO')
	im = ax.imshow(mean(mscs[1].c0s,axis=0).T,vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',
		rotation=270)
	#---plot average structure
	vmax = max([mean(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	vmin = min([mean(msets[i].surf,axis=0).min()/msets[i].lenscale for i in range(2)])
	extrem = max(abs(vmax),abs(vmin))
	vmax,vmin = extrem,-extrem
	if extremz_checkplot_global != None: vmin,vmax = extremz_checkplot_global
	ax = plt.subplot(gs[1,0])
	axeslist.append(ax)
	im = ax.imshow(mean(msets[0].surf,axis=0).T/msets[0].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	ax = plt.subplot(gs[1,1])
	axeslist.append(ax)
	im = ax.imshow(mean(msets[1].surf,axis=0).T/msets[1].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle (\mathrm{nm})$',rotation=270)
	#---plot standard deviations
	extrem = max([std(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
	if extrems_checkplot_global != None: extrem = extrems_checkplot_global
	ax = plt.subplot(gs[2,0])
	axeslist.append(ax)
	ax.imshow((std(msets[0].surf,axis=0)/msets[0].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.jet,
		interpolation='nearest',origin='lower')
	ax = plt.subplot(gs[2,1])
	axeslist.append(ax)
	im = ax.imshow((std(msets[1].surf,axis=0)/msets[1].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.jet,
		interpolation='nearest',origin='lower')
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(r'$\left\langle \left(z-\overline{z}\right)^{2} \right\rangle (\mathrm{{nm}^2})$',
		rotation=270)
	#---compare 2D energy spectra
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		mset = msets[m]
		axr2 = plt.subplot(gs[3,m])
		cm,cn = [int(i/2) for i in shape(mscs[m].tsum2d)]
		wid = 3
		dat = mscs[m].tsum2d[cm-wid:cm+wid+1,cn-wid:cn+wid+1]
		dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
		im = plotter2d(axr2,mset,dat=dat,
			cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
			ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=10,label_style='q',lims=[0.1,10],
			tickskip=int(round(mset.griddims[0]/6,-1)))
	axins = inset_axes(axr2,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=axr2.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(
		r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
		fontsize=fsaxlabel,rotation=270)
	#---compare height-curvature correlation in 2D
	for m in [analysis_names.index(aname) for aname	in plot_reord]:
		mset = msets[m]
		axr3 = plt.subplot(gs[4,m])
		cm,cn = [int(i/2) for i in shape(mscs[m].t2d[1])]
		wid = 3
		dat = mscs[m].t2d[1][cm-wid:cm+wid+1,cn-wid:cn+wid+1]
		dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
		im = plotter2d(axr3,mset,dat=dat,
			cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
			ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=10,label_style='q',lims=[10**-5,10**-2],
			tickskip=int(round(mset.griddims[0]/6,-1)))
	axins = inset_axes(axr3,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=axr3.transAxes,
		borderpad=0)
	axeslist.append(axins)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	axins.set_ylabel(
		r'$\left\langle C_{0,q} h_{-q} \right\rangle $',
		fontsize=fsaxlabel,rotation=270)
	axeslist.append(axins)
	for ax in axeslist:
		plt.setp(ax.get_yticklabels(),fontsize=10)
		plt.setp(ax.get_xticklabels(),fontsize=10)	
	#---spectrum plot
	spectrum_summary(fig=fig,gs=gs2,
		titletext=r'$\mathrm{C_{0,hypo}}='+('{0:.3f}'.format(c0ask))+'a_0^{-1}'+\
		'='+('{0:.3f}'.format(c0ask*msets[1].lenscale))+'\:\mathrm{({nm}^{-1})}$')
	#---save
	plt.savefig(pickles+'fig-bilayer-couple-'+bigname+'.png',bbox_inches='tight',dpi=300)
	if showplots: plt.show()
	plt.close(fig)
	
	'''
	Note: make cool gifs with the following commands on light.site
	gifmake ~/fig-bilayer-couple-compare-v614-120000-220000-200-v2005.gif $(filefield fig-bilayer-couple-compare-v614-120000-220000-200-v2005-C_0-\*png C_0)
	gifmake ~/fig-bilayer-couple-v614-120000-220000-200-v2005.gif $(filefield fig-bilayer-couple-v614-120000-220000-200-v2005-C_0-\*png C_0)		
	'''


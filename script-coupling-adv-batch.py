#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')

import psycopg2
import psycopg2.extras

#---NOTES
#-------------------------------------------------------------------------------------------------------------

"""

check database for simulations
check database for structure pickles
	if no structure pickles, make them
check database for comparisons
	for each missing comparison, make the comparison and then save it
retrieve the comparisons

prepare list of comparisons that you want, from mesosims table
check for the comparisons

"""

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
	'v616-210000-310000-200': {
		'simtype':'md',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'detail_name':r'$\mathrm{{ENTH}\ensuremath{\times}8}$',
		'testname':'v614',
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
	
#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---selection criteria
select_criteria_meso = {'callsign':'v2015'}
cgmd_avail = [
	'v550-300000-400000-200',
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	]

batch_override = True
routine = ['calc','masterplot']
showplots = False

#---naming the batch calculation
mesoname = '-'.join([k for l in [[i,select_criteria_meso[i]] 
	for i in select_criteria_meso.keys()] for k in l])
batchname = \
	'-'.join(cgmd_avail)+'-'+\
	'-'.join([i+'-'+select_criteria_meso[i] for i in select_criteria_meso.keys()])
	
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
			for t in ([0,1,2,3] if (analysis_descriptors[k])['hascurv'] else [0])]) 
			for k in [analysis_names[1],analysis_names[1]]]),
		2*max([max([max([i for i in mscs[analysis_names.index(k)].t1d[t][:,1] if i != 0.]) 
			for t in ([0,1,2,3] 
			if (analysis_descriptors[k])['hascurv'] else [0])]) 
			for k in [analysis_names[1],analysis_names[1]]])))
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
	
def compute_undulation_residuals(simnames,thresh=0.3,specnum=0,view_resids=False):
	'''Compute residuals between undulations for two simulations.'''
	subj = [[],[]]
	xdat = [[],[]]
	ydat = [[],[]]
	for s in range(len(simnames)):
		simname = simnames[s]
		#---use simulation indices to check CGMD vs MESO
		if (analysis_descriptors[simname])['simtype'] == 'md':
			#---the master_spectrum_dict has lots of redundancy
			#---? EXPLAIN
			simname = [i for i in master_spectrum_dict.keys() if re.search(simname,i)][0]
			subj[s] = master_spectrum_dict[simname]
			xdat[s] = collapse_spectrum(subj[s]['cgmd_qs'],subj[s]['cgmd_qs'])
			ydat[s] = collapse_spectrum(subj[s]['cgmd_qs'],subj[s]['cgmd_t2d'][specnum])	
		#---perform mesoscale simulation lookup
		elif (analysis_descriptors[simname])['simtype'] == 'meso':
			subj[s] = master_spectrum_dict[simname]
			xdat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_qs'])
			ydat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_t2d'][specnum])
		else:
			raise Exception('except: no simulation found')
	dat0log = array([i for i in array([xdat[0],log10(ydat[0])]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat[1],log10(ydat[1])]).T if i[0] < thresh])
	resid = mean(abs(dat1log-dat0log)**2)
	tmp = dat1log[:,1]+mean(dat0log[:,1])-mean(dat1log[:,1])
	resid_shift = mean((dat0log-transpose([dat1log[:,0],tmp]))**2)
	if view_resids:
		ax = plt.subplot(111)
		plt.plot(dat0log[:,0],dat0log[:,1],'ro-')
		plt.plot(dat1log[:,0],dat1log[:,1],'bo-')
		plt.show()
	return resid,resid_shift
	
def plotter_undulation_residuals(thresholds,comp,comp_cgmd,toptitle,a0,filename_descriptor):
	'''Plot a summary of the undulation residuals across systems.'''
	fig = plt.figure(figsize=((4,3*len(thresholds)) if len(thresholds)>1 else (8,8)))
	#---?need to generate meso_c0s or ordering from meso_keys here
	sortvals_names = {
		'C_0':r'$\mathrm{C_0}=$',
		'kappa':r'$\mathrm{\kappa}=$',
		}
	meso_names = [' '.join([sortvals_names[sortvals[j]]+str(i[1][j]) 
		for j in range(len(sortvals))]) for i in reord]
	gs = gridspec.GridSpec(len(thresholds),2,wspace=0.1,hspace=0.1,
		width_ratios=[len(meso_keys),len(cgmd_list)])
	cmap = mpl.cm.jet
	#---plots loop over rows where each row has the comparison for a different wavevector cutoff
	for thresh in thresholds:
		residual_comparisons = comp[thresholds.index(thresh)]
		residual_comparisons_cgmd = comp_cgmd[thresholds.index(thresh)]
		axl = plt.subplot(gs[thresholds.index(thresh),0])
		axr = fig.add_subplot(gs[thresholds.index(thresh),1])
		max_resid = max([residual_comparisons.max(),residual_comparisons_cgmd.max()])
		im = axl.imshow(residual_comparisons,interpolation='nearest',cmap=cmap,
			vmax=max_resid,vmin=0)
		axl.set_xticks(range(len(meso_keys)))
		axl.set_yticks(range(len(meso_keys)))
		axl.set_xticklabels(meso_names)
		axl.set_yticklabels(meso_names)
		for label in im.axes.xaxis.get_ticklabels():
			label.set_rotation(90)
		axl.set_xlim((0,len(meso_keys)))
		im = axr.imshow(residual_comparisons_cgmd,interpolation='nearest',cmap=cmap,
			vmax=max_resid,vmin=0)
		#---plot lowest residual rankings
		tmp = array(comp_cgmd)[0].T
		for k in range(4):
			for i in range(len(tmp)):
				pos = argsort(tmp[i])[k]
				axr.text(i,pos,str(k),fontsize=16,horizontalalignment='center',va='center',weight='bold',
					color='w')
		#---settings
		axr.set_xticks(range(len(cgmd_list)))
		axr.set_xticklabels([(analysis_descriptors[name])['label'] for name in cgmd_list])
		for label in im.axes.xaxis.get_ticklabels():
			label.set_rotation(90)
		plt.setp(axr.get_yticklabels(), visible=False)
		axins = inset_axes(axr,width=1./len(cgmd_list)/2.,height="100%",loc=3,
			bbox_to_anchor=(1.2,0.0,1.,1.),
			bbox_transform=axr.transAxes,
			borderpad=0)
		boundaries = linspace(0,max_resid,21)
		ticks = [float('{:.3f}'.format(i)) for i in boundaries]
		norm = mpl.colors.BoundaryNorm(ticks,cmap.N)
		cbar = mpl.colorbar.ColorbarBase(axins,cmap=cmap,ticks=ticks,boundaries=boundaries,
			norm=norm,spacing='proportional',format='%03f')
		cbar.set_ticklabels([str(i) for i in ticks])
		axl.set_xlim((-0.5,len(meso_keys)-1+0.5))
		axr.set_xlim((-0.5,len(cgmd_list)-1+0.5))
		axl.set_title(toptitle,fontsize=fsaxtitle)
		axl.set_ylabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
		axl.set_xlabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
	#---save
	plt.savefig(pickles+'fig-bilayer-couple-meta-'+\
		filename_descriptor+'-'+batchname+'.png',
		bbox_inches='tight',dpi=300)
	if showplots: plt.show()
	plt.close(fig)

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

#---see if the coupling sweep has already been done
master_spectrum_dict = unpickle(pickles+'pkl.bilayer-coupling-sweep.'+batchname+'.pkl')

#---perform the sweep
if master_spectrum_dict == None:
	#---empty dictionary for cross-simulation comparisons
	master_spectrum_dict = {}
	#---loop over all comparisons
	for batch_cgmd in cgmd_avail:
		for batch_meso in [i for i in analysis_descriptors.keys() if i not in cgmd_avail]:
			print '\n'
			match_scales = analysis_names = plot_reord = [batch_cgmd,batch_meso]
			status('status: comparing '+batch_cgmd+' '+batch_meso)
			bigname = '-'.join(analysis_names)
			execfile('script-coupling-adv.py')
			del msets,mscs,collect_c0s
	if master_spectrum_dict != {}:
		pickledump(master_spectrum_dict,'pkl.bilayer-coupling-sweep.'+batchname+'.pkl',directory=pickles)

#---residual sweep parameters
cgmd_list = ['v550-300000-400000-200','v614-120000-220000-200','v616-210000-310000-200']
qmagfilter = (analysis_descriptors[cgmd_avail[0]])['qmagfilter']
#---note: using a threshold of 0.4 does a worse job of matching the kappas, which I can glean by eye
thresholds = [qmagfilter[1],0.4,0.2,0.25][-1:]
thresholds = [0.4,0.3,0.25,0.2][:1]
shift_curves = False

meso_keys = [i for i in analysis_descriptors.keys() 
	if (analysis_descriptors[i]['simtype']=='meso')]
a0 = 1./(master_spectrum_dict[meso_keys[0]])['lenscale']
#---generate a list of curvatures in the mesoscale simulations
if 0 or 'print_new_c0_vals' in routine:
	lister = [0.002, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02, 0.022, 
		0.025, 0.028, 0.03, 0.035, 0.04, 0.045, 0.05]
	lister = [0.002, 0.004, 0.006, 0.008, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 
		0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.032, 0.034, 
		0.036, 0.038, 0.04, 0.045, 0.05]
	lister = [0.01,0.018,0.024,0.030]
	print '( \''+'\' \''.join(['{0:.4f}'.format(round(i*a0,4)) for i in lister])+'\' )'

#---loop over comparisons
if 'comp' not in globals():
	#---perform the fit for both undulations and height-curvature correlations
	specnums = [0]
	fnames = ['undulations','curvature_undulations']	
	toptitles = ['Undulation residuals','Curvature-undulation residuals']
	#---loop
	comp = [[[] for i in thresholds] for sn in specnums]
	comp_cgmd = [[[] for i in thresholds] for sn in specnums]
	#---loop over spectrum comparisons
	for sn in range(len(specnums)):	
		#---loop over desired thresholds for doing the spectrum comparisons
		for thresh in thresholds:
			status('status: threshold = '+str(thresh))
			#---prepare the list of mesoscale simulations for comparison and sort them here
			#---this assumes that the items in analysis_descriptors equal those in master_spectrum_dict
			meso_keys = [i for i in analysis_descriptors.keys() 
				if (analysis_descriptors[i]['simtype']=='meso')] # and analysis_descriptors[i]['C_0'] > 0
			cgmd_keys = [i for i in analysis_descriptors.keys() if analysis_descriptors[i]['simtype']=='md']
			#---sorting the keys
			sortvals = ['kappa','C_0']
			valpairs = [[(analysis_descriptors[i])[j] for j in sortvals] for i in meso_keys]
			reord = sorted(enumerate(valpairs),key=operator.itemgetter(1))
			meso_keys = [meso_keys[i[0]] for i in reord]
			#---compare mesoscale to itself			
			residual_comparisons = zeros((len(meso_keys),len(meso_keys)))
			for i in meso_keys:
				status('status: undulation comparison '+i,i=meso_keys.index(i),looplen=len(meso_keys))
				for j in meso_keys:
					resid,resid_shift = compute_undulation_residuals([i,j],
						thresh=thresh,specnum=sn)
					residual_comparisons[meso_keys.index(i),meso_keys.index(j)] = \
						(resid if not shift_curves else resid_shift)
			#---compare CGMD to mesoscale
			residual_comparisons_cgmd = zeros((len(meso_keys),len(cgmd_list)))
			for i in meso_keys:
				status('status: undulation comparison '+i,i=meso_keys.index(i),looplen=len(meso_keys))
				for name in cgmd_list:
					resid,resid_shift = compute_undulation_residuals([i,name],
						thresh=thresh,specnum=sn)
					residual_comparisons_cgmd[meso_keys.index(i),cgmd_list.index(name)] = \
						(resid if not shift_curves else resid_shift)
			comp[sn][thresholds.index(thresh)] = residual_comparisons
			comp_cgmd[sn][thresholds.index(thresh)] = residual_comparisons_cgmd
		a0 = 1./(master_spectrum_dict[meso_keys[0]])['lenscale']
		plotter_undulation_residuals(thresholds,comp[sn],comp_cgmd[sn],toptitles[sn],
			a0,fnames[sn]+('shifted-' if shift_curves else ''))



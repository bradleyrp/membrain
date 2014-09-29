#!/usr/bin/python

from membrainrunner import *
from numpy import *

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---note that these functions were pulled directly from script-coupling-adv-batch.py
#---note that we may wish to streamline this matching procedure from the script-coupling-adv-batch.py version

def fftwrap_special(dat,redundant=1,shift=True):
	'''This function wraps the standard discrete FFT for a system with possible-redundant rows.'''
	trim = -1 if redundant == 1 else None 
	if shift: return fft.fftshift(fft.fft2(array(dat)[:trim,:trim]))
	else: return fft.fft2(array(dat)[:trim,:trim])

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

#---CLASS
#-------------------------------------------------------------------------------------------------------------

class ModeCouple:
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
	def calculate_mode_coupling(self,mset,c0s,**kwargs):
		'''Moved the entire mode-coupling calculation to a single function so it can be repeated easily.'''
		
		#---added flat for redundant trim here, so the original method requires nort = False
		nort = False
		doshift = False
		
		#---compute dimensions
		mset.calculate_undulations(removeavg=kwargs['removeavg'],fitlims=kwargs['fitlims'],
			forcekappa=kwargs['forcekappa'],qmagfilter=kwargs['qmagfilter'])
		grid = mset.griddims
		if nort: m,n = grid[0]-1,grid[1]-1
		else: m,n = grid
		#---if no curvature field, use zeros
		if c0s == []: c0s = zeros((m,n))
		#---units of hqs are unspecified here but then set by autocorr
		print 'status: calculating hqs'
		#hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
		hqs = []
		#---added shift=False to fix ModeSum 2014.09.18
		for i in range(len(mset.surf)): hqs.append(fftwrap_special(mset.surf[i],
			shift=doshift,
			redundant=nort,
			)/double(m*n))
		#del mset.surf
		#---units of incoming c0s must be in the units specified by mset.lenscale
		#---? units need checked
		print 'status: calculating cqs'
		cqs = [fftwrap_special(mset.lenscale*array(c0s[i]))/double(m*n) for i in range(len(c0s))]
		#print 2
		Lx,Ly = mean(mset.vecs,axis=0)[0:2]
		#print 3
		cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
		#print 4
		qmags = mset.lenscale*array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
			for j in range(0,n)] for i in range(0,m)])
		#print 5
		hqsa = autocorr(hqs,direct=0,lenscale=mset.lenscale)
		#print 6
		hqsb = autocorr(hqs,direct=1,lenscale=mset.lenscale)
		#print 7
		center = [cm,cn]
		#print 8
		cqsa = autocorr(cqs,direct=0,lenscale=1.0)
		#print 9
		cqsb = autocorr(cqs,direct=1,lenscale=1.0)
		#print 10
		qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
		#print 11
		mt,nt = shape(hqsa[0])
		#print 12
		area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
		#print 13
		scalefac = mset.undulate_kappa*area
		#print 14
		#---added raw storage of hqs here
		self.hqs = array(hqs)
		self.cqs = array(cqs)
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
		#print 16
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

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------
		
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
		
def construct_hypofield_deprecated(c0ask):
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
	
def construct_hypofield(fieldspecs,mset=None):
	'''
	Construct a hypothetical curvature field.
	Note that this function feeds the field specification to the gauss2d function in the correct order.
	'''
	#---unpack the box size and lenscale (relative to nm) of the mset object
	if mset != None:
		lenscale = mset.lenscale
		vecs = mean(mset.vecs,axis=0)
		m,n = mset.griddims
	else: raise Exception('except: need to set lenscales in construct_hypofield')
	if 'center_x' not in fieldspecs.keys(): center_x = 0.5
	if 'center_y' not in fieldspecs.keys(): center_y = 0.5
	#---order the specs for gauss2d
	params = [0,
		fieldspecs['C_0'],
		vecs[0]*(1+center_x)/lenscale,
		vecs[1]*(1+center_y)/lenscale,
		fieldspecs['sigma_a'],
		fieldspecs['sigma_b'],
		0]
	#---getgrid is xyz points in nm
	getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
		for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
	c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
		for j in range(3*n)] for i in range(3*m)])
	c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
		for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
	return c0hypo
	
def construct_kappafield(fieldspecs,mset=None):
	'''
	Generate a hypothetical kappa (bending rigidity) field according to a dictionary description.
	'''
	#---unpack the box size and lenscale (relative to nm) of the mset object
	if mset != None:
		lenscale = mset.lenscale
		vecs = mean(mset.vecs,axis=0)
		m,n = mset.griddims
	else: raise Exception('except: need to set lenscales')
	if 'center_x' not in fieldspecs.keys(): center_x = 0.5
	if 'center_y' not in fieldspecs.keys(): center_y = 0.5
	cxr,cyr = vecs[0]*(1+center_x)/lenscale,vecs[1]*(1+center_y)/lenscale

	if fieldspecs['type'] == 'disc':
		#---unpack
		fore = fieldspecs['fore']
		back = fieldspecs['back']
		radius = fieldspecs['radius']
		#---function for the disc
		ic = lambda x,y: fore if (x-cxr)**2+(y-cyr)**2 < radius else back
		#---getgrid is xyz points in nm
		getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
			for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
		c0hypo_nopbc = array([[ic(getgrid[i,j,0],getgrid[i,j,1]) 
			for j in range(3*n)] for i in range(3*m)])
		c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
			for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
		return c0hypo
	else: raise Exception('except: cannot make that field')


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
		
		

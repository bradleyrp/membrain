#!/usr/bin/python -i

from ModeCouple import *
execfile('../locations/header.py')

#---? clean up imports and sync with ModeCouple and membrainrunner
import matplotlib.pyplot as plt

def kappafield(fieldspecs,mset=None):
	'''
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
	
	def incirc(kappa_params,x,y):
		if (x-params[2])**2+(y-params[3])**2 < kappa_params[0]: return kappa_params[1]
		else: return kappa_params[0]
	
	#---getgrid is xyz points in nm
	getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
		for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
	c0hypo_nopbc = array([[incirc([24,30,10],getgrid[i,j,0],getgrid[i,j,1]) 
		for j in range(3*n)] for i in range(3*m)])
	c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
		for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
	return c0hypo

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

callsign = [
	'v550-300000-400000-200',
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	][0]
	
#---hypothesis
hypothesis = {
	'kappa':24,
	'sigma':0.0,
	'curvature':{
		'type':'dimple',
		'C_0':0.001,
		'sigma_a':10,
		'sigma_b':10,
		}
	}

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
mset = unpickle(pickles+analysis_descriptors[callsign]['locate'])

#---make the curvature field
hfield = construct_hypofield(hypothesis['curvature'],mset=mset)

#---kappa field
kf = kappafield(hypothesis['curvature'],mset=mset)
#plt.imshow(array(kf).T,origin='lower',interpolation='nearest');plt.show()

lenscale = mset.lenscale
vecs = mean(mset.vecs,axis=0)
m,n = mset.griddims

kqs = fftwrap(mset.lenscale*array(kf))/double(m*n)
kqsa = autocorr([kqs],direct=0,lenscale=mset.lenscale)
m2,n2 = shape(kqs)
#plt.imshow(abs(kqsa)[0],origin='lower',interpolation='nearest');plt.show()
#plt.imshow(imag(kqsa[0]).T,interpolation='nearest',origin='lower');plt.show()

if 1:

	#---compute mode coupling
	msc = ModeCouple()
	msc.calculate_mode_coupling(mset,
		[hfield for i in range(len(mset.surf))],
		**analysis_descriptors[callsign])

	#---compute the hypothetical spectrum according to our theory
	qmagst = msc.qmagst
	area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
		for i in mset.surf_index])/mset.lenscale**2)
	scalefac = hypothesis['kappa']*area
	tsum2d = scalefac*(msc.t2d[0]*qmagst**4-\
		msc.t2d[1]*qmagst**2-msc.t2d[2]*qmagst**2+msc.t2d[3])+\
		1*qmagst**2*msc.t2d[0]*hypothesis['sigma']*area*mset.lenscale**2
	xdat = collapse_spectrum(msc.qmagst,msc.qmagst)
	ydat = collapse_spectrum(msc.qmagst,tsum2d)

if 1:
	cm,cn = [int(round(i/2.-1)) for i in shape(msc.t2d[0])]
	Lx,Ly = mean(mset.vecs,axis=0)[0:2]
	qmags = mset.lenscale*array([[ [(i-cm)/((Lx)/1.)*2*pi,(j-cn)/((Ly)/1.)*2*pi] for j in range(0,n)] for i in range(0,m)])
	


if 0:
	ax = plt.subplot(111)
	ax.plot(xdat,ydat,'-',lw=1.5,c='g',label='dimple')
	ax.set_xlim((0.01,10))
	ax.set_ylim((0.1,10))
	ax.set_aspect(3./2.)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
	ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
	ax.grid(True)
	plt.show()


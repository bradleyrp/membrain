#!/usr/bin/python

import sys,os
from membrainrunner import *
from ModeCouple import *
from numpy import *

#---? TEMPORARY
import matplotlib.pyplot as plt
import matplotlib as mpl

#---CLASS
#-------------------------------------------------------------------------------------------------------------

class ModeSum:

	'''
	The ModeSum class tests hypothetical curvature fields, bending rigidity fields, and surface tensions on a 
	simulation where these parameters are unknown.
	'''

	def __init__(self,hypothesis,callsign,**kwargs):
		self.hypothesis = hypothesis
		self.callsign = callsign
		bigkeylist = ['pickles','analysis_descriptors']
		self.sets = dict()
		for key in bigkeylist: self.sets[key] = kwargs[key]
		self.msc = None
		self.mset = None
		self.hfield = None
		self.kfield = None
		self.kqs = None
		self.kqqp = None
		
	def fftwrap(self,dat,redundant=1,shift=False):
		'''
		This function wraps the standard discrete FFT.\n
		Note that this code does not drop the final (redundant) edge.
		'''
		if shift: return fft.fftshift(fft.fft2(array(dat)))
		else: return fft.fft2(array(dat))

	def ifftwrap(self,dat,redundant=1,shift=False):
		'''
		This function wraps the standard discrete IFFT.\n
		Note that this code does not drop the final (redundant) edge.
		'''
		if shift: return fft.fftshift(fft.ifft2(array(dat)))
		else: return fft.ifft2(array(dat))
		
	def transforms(self,c0s):
		'''
		Perform Fourier transforms of height and curvature fields.
		Note that this code replaces the ModeCouple code which previously computed these quantities for the
		case where bending rigidity is a constant and hence the $q=q^{\prime}$, that is, there are no cross 
		terms.
		'''
		self.hqs = []
		for i in range(len(self.mset.surf)): 
			self.hqs.append(self.fftwrap(array(self.mset.surf[i]))/double(product(self.mset.griddims)))
		self.cqs = []
		for fr in c0s:
			self.cqs.append(self.fftwrap(array(fr))/double(product(self.mset.griddims)))
		
	def compute_modesum_term(self,term,pklname=None):
		'''
		Compute an individual height-curvature correlation term, in full generality, from the Helfrich with
		curvature. These terms are identified by name below. We use an outer product to perform 
		'''

		#---load
		if self.mset == None:
			self.mset = unpickle(self.sets['pickles']+\
			self.sets['analysis_descriptors'][self.callsign]['locate'])
		
		#---make the curvature field
		if self.hfield == None:
			self.hfield = construct_hypofield(self.hypothesis['curvature'],mset=self.mset)
		
		#---parameters
		lenscale = self.mset.lenscale
		vecs = mean(self.mset.vecs,axis=0)
		m,n = self.mset.griddims

		#---transform height and curvature
		self.transforms([self.hfield for i in range(len(self.mset.surf))])
		m2,n2 = shape(self.hqs)[1:]
			
		#---perform the coupling computation
		status('status: computing term = '+term)
		st = time.time()
		termdat = zeros((m2*n2,m2*n2))
		if term == 'hqs_hqps':
			for fr in range(len(self.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
				termdat += real(outer(self.hqs[fr],self.hqs[fr]))
			termdat = termdat / float(len(self.hqs))
		elif term == 'hqs_c0qps':
			for fr in range(len(self.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
				termdat += real(outer(self.hqs[fr],self.cqs[fr]))
			termdat = termdat / float(len(self.hqs))
		elif term == 'c0qs_hqps':
			for fr in range(len(self.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
				termdat += real(outer(self.cqs[fr],self.hqs[fr]))
			termdat = termdat / float(len(self.hqs))
		elif term == 'c0qs_c0qps':
			#---assumes static field
			termdat = real(outer(self.cqs[0],self.cqs[0]))
		else: raise Exception('requested an unclear term')
		status('\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)+' minutes')
		status('status: saving term as h5 binary')
		if pklname == None: pklname = 'pkl.modecouple.'+self.callsign+'.term-'+term+'.h5pkl'
		else: pklname = pklname + '.h5pkl'
		binarysave(termdat,pklname,directory=self.sets['pickles'])
		return 1
		
	def summer(self,hqs_hqps=None,hqs_c0qps=None,c0qs_hqps=None,c0qs_c0qps=None):
	
		'''
		Function which completes the mode coupling solution after loading individual terms.\n
		Note that this must be followed by a plotting or saving function to handle the result.
		'''

		#---load
		if self.mset == None:
			self.mset = unpickle(self.sets['pickles']+\
			self.sets['analysis_descriptors'][self.callsign]['locate'])
		
		#---make the curvature field
		if self.hfield == None:
			self.hfield = construct_hypofield(self.hypothesis['curvature'],mset=self.mset)

		#---transform height and curvature
		self.transforms([self.hfield for i in range(len(self.mset.surf))])
		m2,n2 = shape(self.hqs)[1:]

		#---parameters
		m2,n2 = shape(self.hqs)[1:]
		self.m2,self.n2 = m2,n2
		cm,cn = [int(round(i/2.-1)) for i in shape(self.hqs[0])]
		Lx,Ly = mean(self.mset.vecs,axis=0)[0:2]

		#---make the bending rigidity field
		if self.kfield == None:
			self.kfield = construct_kappafield(self.hypothesis['kappa'],mset=self.mset)
			#---transform the kappa field
			self.kqs = self.fftwrap(self.mset.lenscale*array(self.kfield)/double(m2*n2))
		
		#---parameters
		lenscale = self.mset.lenscale
		vecs = mean(self.mset.vecs,axis=0)
		m,n = self.mset.griddims
		
		#---construct the $kappa(q+q^{prime})$ matrix
		qs = array(meshgrid(arange(m2),arange(n2))).T
		us = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
		usshift = us-1*(us>array([cm,cn]))*array([m2,n2])
		vs = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
		vsshift = vs-1*(vs>array([cm,cn]))*array([m2,n2])
		biginds = array(meshgrid(arange(m2*n2),arange(m2*n2))).T
		monster = usshift[biginds[...,0]]+vsshift[biginds[...,1]]
		self.kqqp = tile(self.kqs,(3,3))[monster[...,0]+m2,monster[...,1]+n2]

		#---show the computation steps
		#---? DEBUG
		self.vs = vs
		self.us = us
		self.qs = qs
		self.usshift = usshift
		self.vsshift = vsshift
		self.biginds = biginds
		self.monster = monster		

		#---? WORKING ON QMAGS
		if 0:
			qmags = self.mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi 
				for j in range(0,n2) for i in range(0,m2)])
		qmagshift = sqrt(sum((self.usshift/(array([Lx,Ly])*self.mset.lenscale/pi))**2,axis=1))
		qmags = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(array([qmagshift]).T))
		self.qmags = qmags

		#---ADD DEMONDEX TO THE CODE AS A DEMO
		#---FINAL CHECK THAT THE QMAGS ARE CORRECT UNDER THE UNIFIED NAMING SCHEME !!!
		
		#---construct the full Helfrich with curvature in matrix form
		status('status: constructing matrix')
		self.full_matrix = ((qmags*qmags)*(qmags*qmags)*hqs_hqps-(qmags*qmags)*hqs_c0qps-\
			(qmags*qmags)*c0qs_hqps+c0qs_c0qps)*real(self.kqqp)

		#---iterate to larger matrix sizes as a speed test and then diagonlize the final result
		status('status: solution tests')
		for lim in [100,200,500,1000,2000,None]:
			st = time.time()
			self.fullans = linalg.eig(abs(self.full_matrix)[:lim,:lim])
			if lim != None: del self.fullans
			status('status: size = '+str(lim)+' duration = '+str(1./60*(time.time()-st))+' minutes')
'''
m3,n3=m2-1,n2-1
inds = array(where(ms.kqqp>1)).T[1]
[ms.qs[(i)/n3,(i)%n3] for i in inds]
ms.usshift[inds[0]],ms.usshift[inds[1]]
[ms.qs[i/n2,i%n2] for i in j for j in [[1386,2970][1539,2883]]]
sqrt(sum((ms.usshift/(array([Lx,Ly])*ms.mset.lenscale))**2,axis=1))
'''


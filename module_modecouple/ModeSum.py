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
	Needs a hypothesis dictionary and a callsign which is a key in the analysis_descriptors dictionary.
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
		
		self.fullans = None
		
	def compute_modesum_term(self,term,pklname=None):

		'''
		Computation function for the ModeSum class which describes the curvature coupling.
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

		#---compute mode coupling
		if self.msc == None:
			msc = ModeCouple()
			self.msc = msc
			self.msc.calculate_mode_coupling(self.mset,
				[self.hfield for i in range(len(self.mset.surf))],
				**self.sets['analysis_descriptors'][self.callsign])
		m2,n2 = shape(self.msc.hqs)[1:]

		if 0:
			#---compute the hypothetical spectrum according to our theory
			qmagst = self.msc.qmagst
			area = double(mean([self.mset.vec(i)[0]*self.mset.vec(i)[1] 
				for i in self.mset.surf_index])/self.mset.lenscale**2)
			scalefac = self.hypothesis['kappa']*area
			tsum2d = scalefac*(self.msc.t2d[0]*qmagst**4-\
				self.msc.t2d[1]*qmagst**2-self.msc.t2d[2]*qmagst**2+self.msc.t2d[3])+\
				1*qmagst**2*self.msc.t2d[0]*self.hypothesis['gamma']*area*self.mset.lenscale**2
			xdat = collapse_spectrum(self.msc.qmagst,self.msc.qmagst)
			ydat = collapse_spectrum(self.msc.qmagst,tsum2d)
			cm,cn = [int(round(i/2.-1)) for i in shape(self.msc.t2d[0])]
			Lx,Ly = mean(self.mset.vecs,axis=0)[0:2]
			qmags = self.mset.lenscale*array([[ [(i-cm)/((Lx)/1.)*2*pi,(j-cn)/((Ly)/1.)*2*pi] 
				for j in range(0,n)] for i in range(0,m)])

		status('status: computing term = '+term)
		st = time.time()
		termdat = zeros((m2*n2,m2*n2))
		if term == 'hqs_hqps':
			for fr in range(len(self.msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.msc.hqs))
				termdat += real(outer(self.msc.hqs[fr],self.msc.hqs[fr]))
			termdat = termdat / float(len(self.msc.hqs))
		elif term == 'hqs_h0qps':
			for fr in range(len(self.msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.msc.hqs))
				termdat += real(outer(self.msc.hqs[fr],self.msc.cqs[fr]))
			termdat = termdat / float(len(self.msc.hqs))
		elif term == 'h0qs_hqps':
			for fr in range(len(self.msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(self.msc.hqs))
				termdat += real(outer(self.msc.cqs[fr],self.msc.hqs[fr]))
			termdat = termdat / float(len(self.msc.hqs))
		elif term == 'h0qs_h0qps':
			#---assumes static field
			termdat = real(outer(self.msc.cqs[0],self.msc.cqs[0]))
		else: raise Exception('requested an unclear term')
		status('\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)+' minutes')
		status('status: saving term as h5 binary')
		if pklname == None: pklname = 'pkl.modecouple.'+self.callsign+'.term-'+term+'.h5pkl'
		else: pklname = pklname + '.h5pkl'
		binarysave(termdat,pklname,directory=self.sets['pickles'])
		return 1
		
	def summer(self,hqs_hqps=None,hqs_h0qps=None,h0qs_hqps=None,h0qs_h0qps=None):
	
		'''
		'''

		#---load
		if self.mset == None:
			self.mset = unpickle(self.sets['pickles']+\
			self.sets['analysis_descriptors'][self.callsign]['locate'])
		
		#---make the curvature field
		if self.hfield == None:
			self.hfield = construct_hypofield(self.hypothesis['curvature'],mset=self.mset)

		#---compute mode coupling
		if self.msc == None:
			msc = ModeCouple()
			self.msc = msc
			self.msc.calculate_mode_coupling(self.mset,
				[self.hfield for i in range(len(self.mset.surf))],
				**self.sets['analysis_descriptors'][self.callsign])

		#---parameters
		m2,n2 = shape(self.msc.hqs)[1:]
		self.m2,self.n2 = m2,n2
		cm,cn = [int(round(i/2.-1)) for i in shape(self.msc.t2d[0])]
		Lx,Ly = mean(self.mset.vecs,axis=0)[0:2]
		qmags = self.mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi 
			for j in range(0,n2) for i in range(0,m2)])

		#---make the curvature field
		if self.kfield == None:
			self.kfield = construct_kappafield(self.hypothesis['kappa'],mset=self.mset)
			self.kqs = fftwrap(self.mset.lenscale*array(self.kfield))/double(m2*n2)
		
		#---parameters
		lenscale = self.mset.lenscale
		vecs = mean(self.mset.vecs,axis=0)
		m,n = self.mset.griddims
		
		#---INCORRECT METHOD
		status('status: constructing kappa')
		if 0: 
			self.kqqp = [[self.kqs[(int(u/m2)+int(v/n2))%m2,(u%m2+v%n2)%n2] for v in range(m2*n2)] 
				for u in range(m2*n2)]
		
		#---CORRECT, SLOW METHOD
		if 0:	
			#---unravel to index couples
			inds = array([[[int(u/n2)+int(v/n2)-cm,u%n2+v%n2-cn] for v in range(m2*n2)] 
				for u in range(m2*n2)])
			self.inds = inds

			#---wavevectors in the box
			if 0: valids = ((1*(inds[...,0]>-m2)+1*(inds[...,1]>-n2)+\
				1*(inds[...,0]<m2)+1*(inds[...,1]<n2))==4)
			else: valids = ((1*(inds[...,0]>0)+1*(inds[...,1]>0)+\
				1*(inds[...,0]<m2)+1*(inds[...,1]<n2))==4)
			self.valids = valids

			#---construct kqqp
			self.kqqp = zeros((m2*n2,m2*n2))
			if 0: kqqp[inds[...,0],inds[...,1]] = 1
			if 0: kqqp[where(valids)] = 1

			#---loop
			st = time.time()
			for u,v in array(where(valids)).T:
				print u,v,int(u/n2)+int(v/n2),u%n2+v%n2
				self.kqqp[u,v] = self.kqs[int(u/n2)+int(v/n2)-cm,u%n2+v%n2-cn]
			status('status: duration = '+str(1./60*(time.time()-st))+' minutes')

			plt.imshow(self.kqqp,interpolation='nearest',origin='lower');plt.show()
		
		#---CORRECT, SLOW METHOD 2 WHICH IS ACTUALLY INCORRECT
		if 0:
			#---unravel to index couples
			inds = array([[[int(u/n2)+int(v/n2)-cm,u%n2+v%n2-cn] for v in range(m2*n2)] 
				for u in range(m2*n2)])
			self.inds = inds

			#---wavevectors in the box
			valids = ((1*(inds[...,0]>-m2)+1*(inds[...,1]>-n2)+\
				1*(inds[...,0]<m2)+1*(inds[...,1]<n2))==4)
			self.valids = valids

		#qs = reshape(array(meshgrid(arange(-cm,cm+1),arange(-cn,cn+1))).T,(m2*n2,2))
		qs = reshape(array(meshgrid(arange(-cm,cm+1),arange(-cn,cn+1))).T,(m2*n2,2))
		qps = -1*qs
		#---following is probably slow
		st = time.time()
		combos = array([[[i,j] for j in qs] for i in qps])
		status('status: duration = '+str(1./60*(time.time()-st))+' minutes')
	
		#---sum the wavevector pairs in each combination
		cs = sum(combos,axis=2)
		#---isolate the valid elements that could still exist in our box
		valids = (1*(cs[...,0]>-cm)+1*(cs[...,0]<cm)+1*(cs[...,1]>-cn)+1*(cs[...,1]<cn))==4		
		self.valids = valids

		#---compute kqqp
		self.kqqp = zeros((m2*n2,m2*n2))
		us,vs = where(valids)
		vis = array(floor(us/float(n2))+floor(vs/float(n2))-cm).astype(int)
		vjs = array(us%n2+vs%n2-cn).astype(int)
		self.vis,self.vjs = vis,vjs #---TMP
		bigkqs = tile(self.kqs,(3,3))
		self.kqqp[us,vs] = real(bigkqs[vis+m2,vjs+cn])
		
		plt.imshow(real(self.kqqp).T,interpolation='nearest',origin='upper',norm=mpl.colors.LogNorm());plt.show()
		plt.imshow(real(self.kqqp).T,interpolation='nearest',origin='upper');plt.show()
		plt.imshow(real(array(self.kqs)).T,interpolation='nearest',origin='upper');plt.show()

		status('status: constructing matrix')
		self.full_matrix = ((qmags*qmags)*(qmags*qmags)*hqs_hqps+(qmags*qmags)*hqs_h0qps+\
			(qmags*qmags)*h0qs_hqps+h0qs_h0qps)*real(self.kqqp)

		status('status: solution tests')
		for lim in [100,200,500,1000,2000,None]:
			st = time.time()
			self.fullans = linalg.eig(abs(self.full_matrix)[:lim,:lim])
			if lim != None: del self.fullans
			status('status: size = '+str(lim)+' duration = '+str(1./60*(time.time()-st))+' minutes')


#!/usr/bin/python

import sys,os
from membrainrunner import *
from ModeCouple import *
from numpy import *

#---? TEMPORARY
import matplotlib.pyplot as plt

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
		
		status('status: constructing kappa')
		self.kqqp = [[self.kqs[(int(u/m2)+int(v/n2))%m2,(u%m2+v%n2)%n2] 
			for v in range(m2*n2)] for u in range(m2*n2)]

		status('status: constructing matrix')
		self.full_matrix = ((qmags*qmags)*(qmags*qmags)*hqs_hqps+(qmags*qmags)*hqs_h0qps+\
			(qmags*qmags)*h0qs_hqps+h0qs_h0qps)*real(self.kqqp)

		status('status: solution tests')
		for lim in [100,200,500,1000,2000,None]:
			st = time.time()
			self.fullans = linalg.eig(abs(self.full_matrix)[:lim,:lim])
			if lim != None: del self.fullans
			status('status: size = '+str(lim)+' duration = '+str(1./60*(time.time()-st))+' minutes')
		
		
		
		#---DEV DEV DEV
		
		#plt.imshow(array(ms.kfield).T,interpolation='nearest',origin='lower');plt.show()

		if 0:
			#---extra terms
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
				for j in range(0,n2)] for i in range(0,m2)])
			
		

#---DEV
#-------------------------------------------------------------------------------------------------------------

#---development code from previous compute_modesum
if 0:
	if 0:
		#---? deprecated looping method is slower than outer product
		#---? can use memory more efficiently if you split this up and delete the objects after pickling
		if 0:
			print 'computing ensemble averaged coupling terms'
			st = time.time()
			self.hqs_hqps = zeros((m2*n2,m2*n2))
			self.hqs_h0qps = zeros((m2*n2,m2*n2))
			self.h0qs_hqps = zeros((m2*n2,m2*n2))
			self.h0qs_h0qps = zeros((m2*n2,m2*n2))
			for v in range(int(m2*n2)):
				status('v = '+str(v),start=st,i=v,looplen=m2*n2)
				for u in range(int(m2*n2)):
					#---deprecated slow method
					if 0: hqs_hqps[u][v] = [msc.hqs[t][u1,u2]*msc.hqs[t][v1,v2] for t in range(tlen)]
					v1,v2 = int(v/n2),v%n2
					u1,u2 = int(u/n2),u%n2
					self.hqs_hqps[u][v] = mean(msc.hqs[:,u1,u2]*msc.hqs[:,v1,v2])
					self.hqs_h0qps[u][v] = mean(msc.hqs[:,u1,u2]*msc.cqs[:,v1,v2])
					self.h0qs_hqps[u][v] = mean(msc.cqs[:,u1,u2]*msc.hqs[:,v1,v2])
					self.h0qs_h0qps[u][v] = mean(msc.cqs[:,u1,u2]*msc.cqs[:,v1,v2])
			print 1.*(time.time()-st)/60.
			print 'pickling'
			pickledump(self.hqs_hqps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'hqs_hqps.pkl')
			print 'pickling'
			pickledump(self.h0qs_hqps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'h0qs_hqps.pkl')
			print 'pickling'
			pickledump(self.hqs_h0qps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'hqs_h0qps.pkl')
			print 'pickling'
			pickledump(self.h0qs_h0qps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'h0qs_h0qps.pkl')
		if 0:
			print 'computing term 1'
			st = time.time()
			self.hqs_hqps = zeros((m2*n2,m2*n2))
			for fr in range(len(msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(msc.hqs))
				self.hqs_hqps += abs(outer(msc.hqs[fr],msc.hqs[fr]))
			self.hqs_hqps = self.hqs_hqps / float(len(msc.hqs))
			print '\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)
			print 'pickling'
			#pickledump(self.hqs_hqps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'.hqs_hqps.pkl')
			#del self.hqs_hqps
		if 0:
			print 'computing term 2'
			st = time.time()
			self.hqs_h0qps = zeros((m2*n2,m2*n2))
			for fr in range(len(msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(msc.hqs))
				self.hqs_h0qps += abs(outer(msc.hqs[fr],msc.cqs[fr]))
			self.hqs_h0qps = self.hqs_h0qps/ float(len(msc.hqs))
			print '\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)
			print 'pickling'
			pickledump(self.hqs_h0qps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'.hqs_h0qps.pkl')
			del self.hqs_h0qps
		if 0:
			print 'computing term 3'
			st = time.time()
			self.h0qs_hqps = zeros((m2*n2,m2*n2))
			for fr in range(len(msc.hqs)):
				status('fr = '+str(fr),start=st,i=fr,looplen=len(msc.hqs))
				self.h0qs_hqps += abs(outer(msc.cqs[fr],msc.hqs[fr]))
			self.h0qs_hqps = self.h0qs_hqps / float(len(msc.hqs))
			print '\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)
			print 'pickling'
			pickledump(self.h0qs_hqps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'.h0qs_hqps.pkl')
			del self.h0qs_hqps
		if 0:
			print 'computing term 4'
			st = time.time()
			self.h0qs_h0qps = zeros((m2*n2,m2*n2))
			self.h0qs_h0qps = abs(outer(msc.cqs[0],msc.cqs[0]))
			print '\nstatus: duration = '+'{0:.1f}'.format(1.*(time.time()-st)/60.)
			print 'pickling'
			pickledump(self.h0qs_h0qps,self.sets['pickles']+'pkl.modecouple.'+self.callsign+'.h0qs_h0qps.pkl')
			del self.h0qs_h0qps


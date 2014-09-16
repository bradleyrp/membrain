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
		
	def fftwrap(self,dat,redundant=1,shift=False):
		'''This function wraps the standard discrete FFT for a system with possible-redundant rows.
		
		NOTE THAT I HAVE TURNED OFF THE REDUNDANT PART!!
		
		'''
		trim = -1
		if shift: return fft.fftshift(fft.fft2(array(dat)[:trim,:trim]))
		else: return fft.fft2(array(dat)[:trim,:trim])

	def ifftwrap(self,dat,redundant=1,shift=False):
		'''
		This function wraps the standard discrete IFFT for a system with possible-redundant rows.

		NOTE THAT I HAVE TURNED OFF THE REDUNDANT PART!!
		
		'''
		trim = -1
		if shift: return fft.fftshift(fft.ifft2(array(dat)[:trim,:trim]))
		else: return fft.ifft2(array(dat)[:trim,:trim])
		
	def cmc_replace(self,c0s):
		'''
		Replacement function for ModeCouple ...
		NEEDS NOTES HERE
		'''
		self.hqs = []
		for i in range(len(self.mset.surf)): 
			self.hqs.append(self.fftwrap(array(self.mset.surf[i]))/double(product(self.mset.griddims)))
		self.cqs = []
		for fr in c0s:
			self.cqs.append(self.fftwrap(array(fr))/double(product(self.mset.griddims)))
		
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
		if 0:
			#---deprecated method from bloated ModeCouple code
			if self.msc == None:
				msc = ModeCouple()
				self.msc = msc
				self.msc.calculate_mode_coupling(self.mset,
					[self.hfield for i in range(len(self.mset.surf))],
					**self.sets['analysis_descriptors'][self.callsign])
			m2,n2 = shape(self.msc.hqs)[1:]
		else:
			self.cmc_replace([self.hfield for i in range(len(self.mset.surf))])
			m2,n2 = shape(self.hqs)[1:]
			

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
		
		if 1:
			if term == 'hqs_hqps':
			
				cm,cn = [int(round(i/2.-1)) for i in shape(self.hqs[0])]
			
				qs = array(meshgrid(arange(m2),arange(n2))).T
				us = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
				usshift = us-1*(us>array([cm,cn]))*array([m2,n2])
				vs = qs.copy()[arange(m2*n2)/n2,arange(m2*n2)%n2]
				vsshift = vs-1*(vs>array([cm,cn]))*array([m2,n2])
				biginds = array(meshgrid(arange(m2*n2),arange(m2*n2))).T
				monster = usshift[biginds[...,0]]+vsshift[biginds[...,1]]
			
				for fr in range(len(self.hqs)):
					status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
					termdat += real(tile(self.hqs[fr],(3,3))[monster[...,0]+m2,monster[...,1]+n2])
					
				termdat = termdat / float(len(self.hqs))
		

		#---DEPRECATED
		if 0:
			if term == 'hqs_hqps':
				for fr in range(len(self.hqs)):
					status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
					termdat += real(outer(self.hqs[fr],self.hqs[fr]))
				termdat = termdat / float(len(self.hqs))
			elif term == 'hqs_c0qps':
				for fr in range(len(self.hqs)):
					status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
					termdat += real(outer(self.cqs[fr],self.hqs[fr]))
				termdat = termdat / float(len(self.hqs))
			elif term == 'c0qs_hqps':
				for fr in range(len(self.hqs)):
					status('fr = '+str(fr),start=st,i=fr,looplen=len(self.hqs))
					termdat += real(outer(self.hqs[fr],self.cqs[fr]))
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
		'''

		#---load
		if self.mset == None:
			self.mset = unpickle(self.sets['pickles']+\
			self.sets['analysis_descriptors'][self.callsign]['locate'])
		
		#---make the curvature field
		if self.hfield == None:
			self.hfield = construct_hypofield(self.hypothesis['curvature'],mset=self.mset)

		#---compute mode coupling
		if 0:
			#---deprecated method from bloated ModeCouple code
			if self.msc == None:
				msc = ModeCouple()
				self.msc = msc
				self.msc.calculate_mode_coupling(self.mset,
					[self.hfield for i in range(len(self.mset.surf))],
					**self.sets['analysis_descriptors'][self.callsign])
			m2,n2 = shape(self.msc.hqs)[1:]
		else:
			self.cmc_replace([self.hfield for i in range(len(self.mset.surf))])
			m2,n2 = shape(self.hqs)[1:]

		#---parameters
		m2,n2 = shape(self.hqs)[1:]
		self.m2,self.n2 = m2,n2
		cm,cn = [int(round(i/2.-1)) for i in shape(self.msc.t2d[0])]
		Lx,Ly = mean(self.mset.vecs,axis=0)[0:2]
		qmags = self.mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi 
			for j in range(0,n2) for i in range(0,m2)])

		#---make the curvature field
		if self.kfield == None:
			self.kfield = construct_kappafield(self.hypothesis['kappa'],mset=self.mset)
			
			#---INTERVENING HERE TO MAKE SURE KQS IS CORRECT !!!
			self.kqs = fftwrap(
				self.mset.lenscale*array(self.kfield)/double(m2*n2),
				redundant=True,shift=True)
		
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

		#---WORKING CODE, penultimate
		if 0:

			#---wavevectors and opposites
			qs = reshape(array(meshgrid(arange(-cm,cm+1),arange(-cn,cn+1))).T,(m2*n2,2))
			qps = -1*qs
			#---combinations of wavevectors
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
			#---get valid indices
			us,vs = where(valids)
			#---convert indices from couples to singles
			self.vis = array(floor(us/float(n2))+floor(vs/float(n2))-cm).astype(int)
			self.vjs = array(us%n2+vs%n2-cn).astype(int)
			#---tile the transformed field
			self.bigkqs = tile(self.kqs,(3,3))
			#---rearrange the tiled field by filtering out the indices
			self.kqqp[us,vs] = real(self.bigkqs[self.vis+m2,vjs+cn])
		
			if 1:
				ds = (cs[...,0]**2+cs[...,1]**2)**0.5
				ax = plt.subplot(131)
				ax.imshow(valids,interpolation='nearest',origin='upper')
				ax = plt.subplot(132)
				ax.imshow(ds,interpolation='nearest',origin='upper')
				plt.show()
			if 1:
				plt.imshow(real(array(self.kqs)).T,interpolation='nearest',origin='upper');plt.show()
				plt.imshow(real(self.kqqp).T,interpolation='nearest',origin='upper');plt.show()
				plt.imshow(real(self.kqqp).T,interpolation='nearest',origin='upper',
					norm=mpl.colors.LogNorm());plt.show()

		#---wavevectors and opposites
		self.qs = reshape(array(meshgrid(arange(-cm,cm+1),arange(-cn,cn+1))).T,(m2*n2,2))
		self.qps = -1*self.qs.copy()
		#---combinations of wavevectors
		st = time.time()
		self.combos = array([[[i,j] for j in self.qs] for i in self.qps])
		status('status: duration = '+str(1./60*(time.time()-st))+' minutes')
		#---sum the wavevector pairs in each combination
		self.cs = sum(self.combos,axis=2)
		#---tile the transformed field to perform the lookup
		self.bigkqs = tile(self.kqs,(3,3))
		#---transformed kappa is looked up in the order of combos
		self.kqqp = self.bigkqs[self.cs[...,0],self.cs[...,1]]

		status('status: constructing matrix')
		self.full_matrix = ((qmags*qmags)*(qmags*qmags)*hqs_hqps+(qmags*qmags)*hqs_c0qps+\
			(qmags*qmags)*c0qs_hqps+c0qs_c0qps)*real(self.kqqp)

		status('status: solution tests')
		for lim in [100,200,500,1000,2000,None]:
			st = time.time()
			self.fullans = linalg.eig(abs(self.full_matrix)[:lim,:lim])
			if lim != None: del self.fullans
			status('status: size = '+str(lim)+' duration = '+str(1./60*(time.time()-st))+' minutes')


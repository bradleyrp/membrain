#!/usr/bin/env python

#---LIBRARIES
#-------------------------------------------------------------------------------------------------------------

#---Basic libraries
import os
if os.path.isfile('/etc/pythonstart'):
	execfile('/etc/pythonstart')
import time
import argparse
import string
import importlib
import os
import errno
import sys
import pickle
import collections
import time
import random
import subprocess
import re
import code
import datetime
import operator
try: import h5py
except: pass

#---Set up visualization
#---? this is likely required for Mayavi, so consider moving to plotter.py
os.environ['ETS_TOOLKIT'] = 'qt4'

#---MDAnalysis functions
#---note that 
from MDAnalysis import *

#---Mathematics libraries
from numpy import *
import numpy as np
import scipy.spatial
import scipy.interpolate

#---Custom classes
from membraindata import *

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---execute membrainrunner.py only to get the correct status function
#---check the parent directory in case we are executing in a subfolder
if os.path.isfile('membrainrunner.py'): execfile('membrainrunner.py')
elif os.path.isfile('../membrainrunner.py'): execfile('../membrainrunner.py')
else: raise Exception('except: cannot locate membrainrunner.py')

#---MEMBRANESET CLASS
#-------------------------------------------------------------------------------------------------------------

class MembraneSet:
	'''A class which holds a membrane simulation.'''
	def __init__(self):
		#---Raw data
		self.xyzs = []
		self.xyzs_time = []
		self.tri = []
		self.tri_index = []
		#---Trajectory descriptors
		self.nframes = 0
		self.resolution = ''
		self.vecs = []
		self.vecs_index = []
		self.griddims = []
		self.universe = []
		self.universe_structfile = ''
		self.universe_trajfile = ''
		self.monolayer_residues = []
		self.monolayer_residues_abs = []
		self.resnames = []
		self.resids = []
		self.resids_abs = []
		self.monolayer_by_resid = []
		self.monolayer_by_resid_abs = []
		self.resids_reracker = []
		self.monolayer_rep = ''
		self.selections = []
		self.picklename = ''
		self.time_start = 0.0
		self.time_total = 0.0
		self.time_dt = 0.0
		self.time_list = []
		#---Surface heights on a grid
		self.surf = []
		self.surf_thick = []
		self.surf_index = []
		self.surf_time = []
		self.monolayer1 = []
		self.monolayer2 = []
		self.surf_mean = []
		self.surf_position = []
		self.surf_positioni = []
		#---Undulation data
		self.undulate_raw = []
		self.undulate_hqhq2d = []
		self.undulate_qmag2d = []
		self.undulate_spec1d = []
		self.undulate_peri_raw = []
		self.undulate_peri_hqhq2d = []
		self.undulate_peri_spec1d = []
		self.undulate_qmagfilter = []
		self.undulate_kappa = 0.0
		#---Protein data
		self.protein = []
		self.protein_index = []
		#---Length scale parameters
		self.rounder = None
		self.lenscale = None
		#---List of MembraneData objects
		self.store = []
	def __getstate__(self):
		'''Returns the state for pickling.'''
		result = self.__dict__.copy()
		del result['universe']
		if (self.resolution != 'aamd' and self.resolution != 'cgmd') or (self.xyzs == []):
			del result['xyzs']
		return result

#---Loading functions

	def load_trajectory(self,files,resolution=None,lenscale=None):
		"""Load a molecular dynamics trajectory into the MembraneSet instance.

        This function takes a structure and trajectory file and creates a Universe instance from the
        MDAnalysis library which serves as a member of MembraneSet. The function also populates meta-data 
        inside the MembraneSet instance which keep track of the frames and timining of the trajectory.

        Args:
           files (tuple): A tuple containing the path to the structure (usually a GRO file) and trajectory for
           the target simulation.

        Kwargs:
           resolution (str): Specify either "CGMD" or "AAMD" to denote the "graining" of simulation. This is
           only used in a few functions, for bookkeeping purposes.
           
           lenscale (float): Set the natural lengthscale of the incoming data. Since this function uses the
           MDAnalysis package, the default value is 10, corresponding to 10 Angstroms, so that all 
           downstream units are in nanometers.
           
        You never call this class before calling :func:`public_fn_with_sphinxy_docstring`.

        .. note::

        	An example of intersphinx is this: you **cannot** use :mod:`pickle` on this class.

        """
		
		self.lenscale = lenscale if lenscale != None else 10.
		if resolution != None:
			self.resolution = resolution
		if self.nframes != 0:
			print 'status: clearing trajectory'
			self.vecs = []
			self.vecs_index = []
			self.griddims = []
			self.nframes = 0
		self.universe = Universe(files[0],files[1])
		#---Nb added these flags after fixing transpose error, so they are proxy
		self.universe_structfile = files[0]
		self.universe_trajfile = files[1]
		self.nframes = len(self.universe.universe.trajectory)
		if hasattr(self.universe.trajectory[0],'time'):
			self.time_start = self.universe.trajectory[0].time
		if hasattr(self.universe.trajectory,'totaltime'):
			self.time_total = self.universe.trajectory.totaltime
		if hasattr(self.universe.trajectory,'dt'):
			self.time_dt = self.universe.trajectory.dt
		'''
		if hasattr(self.universe,'trajectory'):
			self.time_list = [self.universe.trajectory[i].time 
				for i in range(len(self.universe.trajectory)) 
				if hasattr(self.universe.trajectory[i],'time')]
		'''
		status('status: the trajectory file has '+str(self.nframes)+' frames')

	def load_points(self,xyzdir,nbase=0,start=None,end=None,
		prefix='conf-',suffix='.xyz',xyzform='',rounder=1.0,lenscale=None,skip=1,regular=False,shifter=None):
		'''Load a trajectory of xyz files.'''
		if self.nframes != 0:
			print 'starus: clearing trajectory'
			self.vecs = []
			self.vecs_index = []
			self.griddims = []
			self.xyzs = []
			self.surf = []
			self.nframes = 0
		self.lenscale = lenscale if lenscale != None else 1.0
		suffix = suffix.replace('\\','')
		dirlistall = os.listdir(xyzdir)
		dirlist = [int(i[len(prefix):-len(suffix)]) for i in dirlistall 
			if i[:len(prefix)] == prefix and i[-len(suffix):] == suffix]
		dirlist.sort()
		if start == None and end == None: 
			whichframes = dirlist[::skip]
		else:
			whichframes = dirlist[dirlist.index(start):dirlist.index(end):skip]
		print 'status: trajectory directory has '+str(len(whichframes))+' frames'
		dataset = []
		if regular == True:
			self.vecs = [[int((nbase-1)*lenscale*rounder),int((nbase-1)*lenscale*rounder)] \
				for i in range(len(whichframes))]
			self.vecs_index = range(len(whichframes))
			self.griddims = [int(nbase),int(nbase)]
			for index in whichframes:
				print 'status: loading '+str(index)
				fp = open(xyzdir+'/'+prefix+str(index)+suffix,'r')
				frame_temp = []
				for line in fp:
					frame_temp.append([float(i) for i in line.split()])
				self.xyzs.append(array(frame_temp))
				self.surf.append(self.rezipgrid(self.xyzs[-1],diff=0))
				self.nframes += 1
		elif regular == False:
			if xyzform == 'rect':
				for index in whichframes:
					fp = open(xyzdir+'/'+prefix+str(index)+suffix,'r')
					frame_temp = []
					for line in fp:
						frame_temp.append([float(i) for i in line.split()])
					self.xyzs.append(array(frame_temp))
					self.nframes += 1
					self.rounder = rounder
				scalefac = 1.3
				self.vecs = [[scalefac*sqrt(3)/2*int(nbase),scalefac*int(nbase)] for index in whichframes]
				self.vecs_index = range(len(whichframes))
				self.griddims = [int(round(self.vecs[0][0]/rounder)),int(round(self.vecs[0][1]/rounder))]
			elif xyzform == 'square':	
				shift = [-shifter,-shifter,0.]
				for index in whichframes:
					fp = open(xyzdir+'/'+prefix+str(index)+suffix,'r')
					frame_temp = []
					for line in fp:
						frame_temp.append([float(line.split()[i])+shift[i] for i in range(len(line.split()))])
					self.xyzs.append(array(frame_temp))
					self.nframes += 1
					self.rounder = rounder
				self.vecs = [[int((nbase-1)*rounder*lenscale),int((nbase-1)*rounder*lenscale)] \
					for index in whichframes]
				self.vecs_index = range(len(whichframes))
				self.griddims = [int(nbase),int(nbase)]
			elif xyzform == 'square2':	
				shift = [-shifter,-shifter,0.]
				for index in whichframes:
					fp = open(xyzdir+'/'+prefix+str(index)+suffix,'r')
					frame_temp = []
					firstline = ''
					for line in fp:
						if firstline == '':
							firstline = line
						else:
							frame_temp.append([float(line.split()[i])+shift[i] for i in 
								range(len(line.split()))])
					self.xyzs.append(array(frame_temp))
					self.nframes += 1
					self.rounder = rounder
				self.vecs = [[int((nbase-1)*rounder*lenscale),int((nbase-1)*rounder*lenscale)] \
					for index in whichframes]
				self.vecs_index = range(len(whichframes))
				self.griddims = [int(nbase),int(nbase)]
		
	def load_points_vtu(self,vtudir,extra_props=None,nbase=0,start=None,end=None,
		prefix='equib-conf-',suffix='.vtu',xyzform='',rounder=1.0,lenscale=None,skip=1):
		'''Load points and extra data from paraview-style VTU files.'''
		import xml.etree.ElementTree as ET
		if self.nframes != 0:
			print 'status: clearing trajectory'
			self.vecs = []
			self.vecs_index = []
			self.griddims = []
			self.xyzs = []
			self.surf = []
			self.nframes = 0
		self.lenscale = lenscale if lenscale != None else 1.0
		suffix = suffix.replace('\\','')
		dirlistall = os.listdir(vtudir)
		dirlist = [int(i[len(prefix):-len(suffix)]) for i in dirlistall 
			if i[:len(prefix)] == prefix and i[-len(suffix):] == suffix]
		dirlist.sort()
		if start == None and end == None: 
			whichframes = dirlist[::skip]
		else:
			whichframes = dirlist[dirlist.index(start):dirlist.index(end):skip]
		print 'status: trajectory directory has '+str(len(whichframes))+' frames'
		extradat = []
		boundary_pts = []
		for index in whichframes:
			tree = ET.parse(vtudir+'/'+prefix+str(index)+suffix)
			root = tree.getroot()
			if index == whichframes[0]:
				extra_item_nums = []
				if extra_props != None:
					if type(extra_props) == str: extra_props = [extra_props]
					for propstring in extra_props:
						extra_item_nums.append([j[1][1] for j in [i.items() 
							for i in root[0].find('Piece').find('PointData').\
							getchildren()]].index(propstring))
				boundsind = [j[1][1] for j in [i.items() for i in root[0].find('Piece').find('PointData').\
					getchildren()]].index('boundary')
				nonbounds = list(where(array(map(float,root[0].find('Piece').find('PointData').\
					getchildren()[boundsind].text.split()))==0.0)[0])
			status('status: reading vtu file number',i=whichframes.index(index),looplen=len(whichframes))
			coord = root[0].find('Piece').find('Points').find('DataArray').text.split()
			coord = map(float,coord)
			n = coord.__len__()/3
			vcoord = (numpy.array(coord).reshape(n,3))
			somedata = []
			for itemno in extra_item_nums:
				somedata.append(map(float,root[0].find('Piece').find('PointData').\
					getchildren()[itemno].text.split()))
			if extra_props != None:
				extradat.append(array(somedata)[:,nonbounds])
			self.xyzs.append(array(vcoord)[nonbounds])
			self.nframes += 1
			self.rounder = rounder
		scalefac = 1.3
		self.vecs = [[scalefac*sqrt(3)/2*int(nbase),scalefac*int(nbase)] for index in whichframes]
		self.vecs_index = range(len(whichframes))
		self.griddims = [int(round(self.vecs[0][0]/rounder)),int(round(self.vecs[0][1]/rounder))]
		return extradat
	
	def gotoframe(self,frameno):
		'''Iterate to another frame, quickly.'''
		if len(self.universe.trajectory) == 1:
			print 'warning: only one frame is available'
		elif frameno == 0:
			frame = self.universe.trajectory[frameno]
		elif self.universe.trajectory.frame-1 < frameno:
			[self.universe.trajectory.next() for i in range(frameno-(self.universe.trajectory.frame-1))]
		elif self.universe.trajectory.frame-1 > frameno:
			frame = self.universe.trajectory[frameno]

#---General identification functions

	def get_points(self,frameno,selection_index=-1):
		'''Shell function that returns coordinates for a pre-defined selection.'''
		if self.universe.trajectory.frame-1 != frameno:
			self.gotoframe(frameno)
		return self.selections[selection_index].coordinates()
		
	def vec(self,frameno):
		'''Shell function that returns box vectors.'''
		if frameno not in self.vecs_index:
			if self.universe.trajectory.frame-1 != frameno:
				self.gotoframe(frameno)
			vec = self.universe.dimensions[0:3]
			self.vecs.append(vec)
			self.vecs_index.append(frameno)
			return vec
		else:
			return self.vecs[self.vecs_index.index(frameno)]

	def identify_monolayers(self,atomdirectors,startframeno=0):
		'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
		status('status: identifying monolayers')
		if 0: status('status: moving to frame '+str(startframeno))
		self.gotoframe(startframeno)
		pointouts = self.universe.selectAtoms(atomdirectors[0])
		pointins = [self.universe.selectAtoms(atomdirectors[j]).coordinates() 
			for j in range(1,len(atomdirectors))]
		whichlayer = [0 if i > 0.0 else 1 for i in pointouts.coordinates()[:,2] - mean(pointins,axis=0)[:,2]]
		monos = []
		#---monos separates the lipids by absolute index into monolayers from index zero
		monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 0])
		monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 1])
		#---monolayer rerack hack if some residue IDs are missing
		#---Nb this may affect the tilter, mesher, identify_residues, and batch_gr functions so beware
		if (max(monos[0]+monos[1])-min(monos[0]+monos[1])) != len(monos[0]+monos[1])-1:
			print 'warning: resorting the monolayer indices because there is a mismatch'
			#---reracker is a sorted list of all of the absolute indices
			reracker = list(sort(monos[0]+monos[1]))
			#---monos_rerack is a copy of reracker in relative indices which is separated into monolayers 
			monos_rerack = [[reracker.index(i) for i in monos[m]] for m in range(2)]
			self.monolayer_residues_abs = monos
			#---resids_reracker is in absolute units
			self.resids_reracker = reracker
			#---monolayer_residues is in relative units
			self.monolayer_residues = monos_rerack
		else: self.monolayer_residues = monos
		if len(monos[0]) != len(monos[1]):
			print 'warning: there is a difference in the number of lipids per monolayer'

	def identify_residues(self,selector):
		'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
		self.monolayer_by_resid = []
		self.resnames = []
		self.resids = []
		#---selector holds the proposed resnames
		for sel in selector:
			selection = self.universe.selectAtoms('resname '+sel)
			if len([i-1 for i in selection.resids()]):
				#---resids holds the absolute residue numbers
				self.resids.append([i-1 for i in selection.resids()])
				self.resnames.append(sel)
		#---when re-racker has been defined by identify_monolayers then provide distinct residue numberings
		if self.resids_reracker != []: 
			self.resids_abs = self.resids
			newresids = [[self.resids_reracker.index(i) for i in j] for j in self.resids]
			self.resids = newresids
			#self.resids = [[flatten(self.monolayer_residues).index(i) for i in j] for j in self.resids]

		else:
			self.resids_abs = self.resids
		#---monolayer_residues is in relative indices
		for monolayer in self.monolayer_residues:
			monolayer_resids = []
			#for resids in residue_ids:
			#---self.resids is in relative indices
			for resids in self.resids:
				monolayer_resids.append(list(set.intersection(set(monolayer),set(resids))))
			#---monolayer_by_resid is in relative indices
			self.monolayer_by_resid.append(monolayer_resids)
		if self.resids_reracker == []:
			self.monolayer_by_resid_abs = self.monolayer_by_resid
		else:
			#---pack absolute residue indices into monolayer_by_resid_abs
			self.monolayer_by_resid_abs = \
				[[[self.resids_reracker[r] for r in restype] 
				for restype in mono] for mono in self.monolayer_by_resid]

	def locate_bilayer(self,frameno,monolayer_rep=None):
		'''General side-of monolayer identifier function. Needs name of the atom.'''
		if frameno not in self.surf_positioni:
			if self.universe.trajectory.frame-1 != frameno:
				self.gotoframe(frameno)
			if monolayer_rep == None: monolayer_rep = self.monolayer_rep
			selection = self.universe.selectAtoms('name '+monolayer_rep)
			position = mean(array(selection.coordinates())[:,2])
			self.surf_position.append(position)
			self.surf_positioni.append(frameno)
			return position
		else:
			return self.surf_position[self.surf_positioni.index(frameno)]
			
#---Shape functions

	def surfacer(self,skip=1,interp='best',lenscale=None):
		'''Interpolate the mesoscale bilayers.'''
		for fr in range(0,len(self.xyzs),skip):
			status('status: interpolating splines',i=fr,looplen=len(range(0,len(self.xyzs),skip)))
			#---Nb currently set for VTU files from RWT simulations. Needs explained.
			vecs = self.vecs[fr]
			grid = self.griddims[0],self.griddims[1]
			vecs = [vecs[i]-vecs[i]/(grid[i]-1) for i in range(2)]
			res0 = self.wrappbc(array(self.xyzs[fr]),vecs=vecs,mode='grow')
			res1 = self.makemesh(res0,vecs,self.griddims,method=interp)
			rezip = self.rezipgrid(res1,diff=1)
			self.surf.append(rezip-mean(rezip))
			self.surf_index.append(fr)
		#---scale box vectors by the correct length scale, deprecated
		if 0:
			self.vecs = [[j*self.lenscale for j in i] for i in self.vecs]
			self.surf_index = range(0,len(self.xyzs),skip)
			
	def surfacer_general(self,surfdata,skip=1,interp='best'):
		'''Interpolate any quantity, as long as you supply it in the order of the xyz positions.'''
		interpdata = []
		for fr in range(0,len(self.xyzs),skip):
			status('status: interpolating splines',i=fr,looplen=len(range(0,len(self.xyzs),skip)))
			#---Nb currently set for VTU files from RWT simulations. Needs explained.
			vecs = self.vecs[fr]
			grid = self.griddims[0],self.griddims[1]
			vecs = [vecs[i]-vecs[i]/(grid[i]-1) for i in range(2)]
			res0 = self.wrappbc(array([[self.xyzs[fr][j][0],self.xyzs[fr][j][1],surfdata[fr][j]] 
				for j in range(len(surfdata[fr]))]),vecs=vecs,mode='grow')
			res1 = self.makemesh(res0,vecs,grid,method=interp)
			rezip = self.rezipgrid(res1,diff=1)
			interpdata.append(rezip)
		return interpdata
			
	def midplaner(self,selector,rounder=4.0,framecount=None,start=None,end=None,skip=None,
		interp='best',protein_selection=None,residues=None,timeslice=None,thick=False,fine=None):
		'''Interpolate the molecular dynamics bilayers.'''
		self.rounder = rounder
		if timeslice != None:
			start = int((float(timeslice[0])-self.time_start)/timeslice[2])
			end = int((float(timeslice[1])-self.time_start+self.time_dt)/self.time_dt)
			skip = int(float(timeslice[2])/self.time_dt)
			skip = 1 if skip < 1 else skip
			print 'status: starting midplaner with [start,end,skip] = ['+\
				str(start)+','+str(end)+','+str(skip)+']'
			if end > len(self.universe.trajectory):
				print 'warning: there are fewer frames than you inferred so dt may be inconsistent'
				end = len(self.universe.trajectory)
				print 'status: starting midplaner with [start,end,skip] = ['+\
					str(start)+','+str(end)+','+str(skip)+']'
		elif framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		self.griddims = [int(round(self.vec(1)[0]/rounder)),int(round(self.vec(1)[1]/rounder))]
		starttime = time.time()
		for k in range(start,end,skip):
			status('status: calculating midplane fr = '+str(k+1)+'/'+str(len(range(start,end,skip))))
			self.calculate_midplane(selector,k,rounder=rounder,interp=interp,residues=residues,
				thick=thick,fine=fine)
			if protein_selection != None:
				self.protein.append(self.universe.selectAtoms(protein_selection).coordinates())
				self.protein_index.append(k)
			
	def mesher(self,selector,framecount=None,start=None,end=None,skip=None,protein_selection=None):
		'''Create a standard mesh from the bilayer surface.'''
		#---Not set up to handle edge cases. Only use for curvature calculation near proteins.
		if framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		for k in range(start,end,skip):
			print 'status: calculating mesh for frame: '+str(k)
			starttime = time.time()
			points = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for i in self.monolayer_residues[0]])
			self.tri.append(scipy.spatial.Delaunay(points[:,0:2]))
			self.xyzs.append(points)
			self.tri_index.append(k)
			if protein_selection != None:
				self.protein.append(self.universe.selectAtoms(protein_selection).coordinates())
				self.protein_index.append(k)
			
	def calculate_midplane(self,selector,frameno,pbcadjust=1,rounder=4.0,interp='best',storexyz=True,
		residues=None,thick=False,fine=None):
		'''Find the midplane of a molecular dynamics bilayer.'''
		#---temporary
		st = time.time()
		lenscale = self.lenscale
		self.gotoframe(frameno)
		#---Use all residues
		#---Note: if you are getting an error below, you probably have redundant residue numbers.
		if residues == None:
			topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for i in self.monolayer_residues[0]])
			botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for i in self.monolayer_residues[1]])
		else:
			topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for r in residues for i in self.monolayer_by_resid[0][self.resnames.index(r)]])
			botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for r in residues for i in self.monolayer_by_resid[1][self.resnames.index(r)]])
		#---new addition for sharp
		#topxyz = topxyz - mean(topxyz)
		#---new addition for sharp
		#botxyz = botxyz - mean(botxyz)
		if storexyz: self.xyzs.append([topxyz,botxyz])
		#---First we wrap PBCs
		topxyzwrap = self.wrappbc(topxyz,vecs=self.vec(frameno),mode='grow')
		botxyzwrap = self.wrappbc(botxyz,vecs=self.vec(frameno),mode='grow')
		#---Triangulate the surface on a regular grid
		topmesh = self.makemesh(topxyzwrap,self.vec(frameno),self.griddims,method=interp,fine=fine)
		botmesh = self.makemesh(botxyzwrap,self.vec(frameno),self.griddims,method=interp,fine=fine)
		#---Take the average surface
		#---Nb this is the location of the infamous transpose errors. Forgot to reverse order of list indices
		topzip = self.rezipgrid(topmesh,frameno=frameno)
		botzip = self.rezipgrid(botmesh,frameno=frameno)

		#plt.imshow(tile(array(topzip).T,(3,3)),interpolation='nearest',origin='lower')
		#plt.show()

		#---Convert from points in 3-space to heights on a grid
		self.monolayer1.append(array(topzip))
		self.monolayer2.append(array(botzip))
		surfz = [[1./2*(topzip[i][j]+botzip[i][j]) for j in range(self.griddims[1])] 
			for i in range(self.griddims[0])]
		if thick:
			surf_thick = [[(topzip[i][j]-botzip[i][j]) for j in range(self.griddims[1])] 
				for i in range(self.griddims[0])]
			self.surf_thick.append(surf_thick)
		self.surf_position.append(mean(surfz))
		surfz = surfz - mean(surfz)
		self.surf.append(surfz)
		self.surf_index.append(frameno)
		self.surf_time.append(self.universe.trajectory[frameno].time)
		#print 1./60*(time.time()-st)

	def triangulator(self,selector,start=None,end=None,skip=None,framecount=None,label=None,tesstype=None):
		'''Triangulate the surface by lipid for the entire trajectory.'''
		if framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		if label == None:
			label = ('cells' if (tesstype == None or tesstype == 'voronoi') else 'triangles')
		result_data = MembraneData(
			('cells' if (tesstype == None or tesstype == 'voronoi') else 'triangles'),label=label)
		for k in range(start,end,skip):
			print 'status: calculating triangulation for frame: '+str(k)
			ans = self.calculate_triangulate(selector,k,tesstype=tesstype)
			result_data.add(ans,[k])
		self.store.append(result_data)

	def calculate_triangulate(self,selector,frameno,tesstype=None):
		'''Produce a clean tesselation.'''
		self.gotoframe(frameno)
		if type(selector) == str:
			selector = [selector for i in range(len(self.resnames))]
		top = []
		for t in range(len(self.resnames)):
			top.extend([mean(self.universe.residues[i].selectAtoms(selector[t]).coordinates(),axis=0) 
				for i in self.monolayer_by_resid[0][t]])
		bot = []
		for t in range(len(self.resnames)):
			bot.extend([mean(self.universe.residues[i].selectAtoms(selector[t]).coordinates(),axis=0) 
				for i in self.monolayer_by_resid[1][t]])
		return [self.calculate_tesselation(array(top),frameno,tesstype=tesstype),
			self.calculate_tesselation(array(bot),frameno,tesstype=tesstype)]

	def calculate_tesselation(self,points,frameno,separate=False,tesstype=None):
		'''Generate the tesselation according to desired type.'''
		if tesstype == None: tesstype = 'voronoi'
		if tesstype == 'voronoi':
			points_pbc = self.wrappbc(points,vecs=self.vec(frameno),mode='grow')
			vmap = scipy.spatial.Voronoi(points_pbc[:,0:2])
			areas = []
			for p in range(len(points)):
				vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
				pairs = zip(vertices, vertices[1:] + vertices[0:1])
				ans = abs(sum(x1 * y2 - y1 * x2 for (x1, y1), (x2, y2) in pairs) / 2)
				areas.append(ans)
			result = [points_pbc,vmap,areas]
		#---Other tesselation methods are deprecated by the Voronoi method above, which is far superior.
		else:
			points_pbc = self.wrappbc(points,vecs=self.vec(frameno),mode='grow')
			pointspbctri = Delaunay(points_pbc[:,0:2])
			pointstri = Delaunay(points[:,0:2])
			simpspbc = [list(i) for i in pointspbctri.simplices]
			simps = [list(i) for i in pointstri.simplices]
			simpsi_g1 = list(where(array([sum([simpspbc[i][j]<400 for j in range(3)]) 
				for i in range(len(simpspbc))])>=2)[0])
			simpsi_g2 = list(where(array([sum([simpspbc[i][j]<400 for j in range(3)]) 
				for i in range(len(simpspbc))])==1)[0])
			simps_g1 = [simpspbc[i] for i in simpsi_g1]
			simps_g2 = [simpspbc[i] for i in simpsi_g2]
			points_g1 = list(set([i for j in simps_g1 for i in j]))
			points_g2 = list(set([i for j in simps_g2 for i in j]))
			#---Tesselate perfectly, irrespective of area considerations (for plotting)
			if tesstype == 'perfect':
				moves = [[0,1,0],[1,0,0],[-1,0,0],[0,-1,0]]
				lost = []
				for tri in simps_g2:
					copies = [[i+move*self.vec(frameno) for i in [points_pbc[i] 
						for i in tri]] for move in moves]
					vt = [[i[0] for i in j] for j in [[list(where(all(points_pbc==i,axis=1))[0]) for i in j] 
						for j in copies] if all([len(k)!=0 for k in j])]
					if not any([set(j) in [set(i) for i in simps_g1] for j in vt]):
						lost.append(vt)
				lost = [i for j in lost for i in j]
				lost = [k for k in lost if sum([[t in points_g1] for t in k])==2]
				lost = [list(j) for j in list(set(tuple([frozenset(i) for i in lost])))]
				found = [lost[l] for l in list(where([i > self.vec(frameno)[0]/2. 
					for i in [mean([points_pbc[i] for i in j],axis=0)[0] for j in lost]])[0])]
				validpts = list(set(points_g1+[i for j in found for i in j]))
				validsimps_unsorted = simps_g1+found
				validsimps = [[validpts.index(j) for j in i] for i in validsimps_unsorted]
				result = [array([points_pbc[i] for i in validpts]),validsimps]
			#---Tesselate for area-calculation purposes
			elif tesstype == 'imperfect':
				validpts = list(set(points_g1+points_g2))
				validsimps_unsorted = simps_g1+simps_g2
				validsimps = [[validpts.index(j) for j in i] for i in validsimps_unsorted]
				result = [array([points_pbc[i] for i in validpts]),validsimps]
			#---General tesselation
			elif tesstype == 'both':
				moves = [[0,1,0],[1,0,0],[-1,0,0],[0,-1,0]]
				lost = []
				for tri in simps_g2:
					copies = [[i+move*self.vec(frameno) for i in [points_pbc[i] 
						for i in tri]] for move in moves]
					vt = [[i[0] for i in j] for j in [[list(where(all(points_pbc==i,axis=1))[0]) for i in j] 
						for j in copies] if all([len(k)!=0 for k in j])]
					if not any([set(j) in [set(i) for i in simps_g1] for j in vt]):
						lost.append(vt)
				lost = [i for j in lost for i in j]
				lost = [k for k in lost if sum([[t in points_g1] for t in k])==2]
				lost = [list(j) for j in list(set(tuple([frozenset(i) for i in lost])))]
				found = [lost[l] for l in list(where([i > self.vec(frameno)[0]/2. 
					for i in [mean([points_pbc[i] 
					for i in j],axis=0)[0] for j in lost]])[0])]
				validpts = list(set(points_g1+points_g2))
				validpts = list(set.union(set(validpts),set([i for j in found for i in j])))
				validsimps_unsorted_perfect = simps_g1+found
				validsimps_perfect = [[validpts.index(j) for j in i] for i in validsimps_unsorted_perfect]
				#---We could just dump both two-ghost and one-ghost triangles as follows:
				#---"validsimps_unsorted_g2 = simps_g1+simps_g2" but this is inefficient and hard to parse. 
				#---Check no redundancies in outputted triangles later on.
				#---Also check whether the Voronoi method above uses midpoints or centroids.
				validsimps_unsorted_g2 = [i for i in simps_g2 if i not in validsimps_perfect]
				validsimps_g2 = [[validpts.index(j) for j in i] for i in validsimps_unsorted_g2]
				#---We return the perfect tesselation triangles plus the two-ghost triangles needed for area 
				#---calculations followed by the length of the simplex list necessary to perfectly tesselate
				result = [array([points_pbc[i] for i in validpts]),
					validsimps_perfect+validsimps_g2,
					len(validsimps_perfect)]
		return result

	def calculate_average_surface(self):
		'''Calculates the average surface for constant area surfaces.'''
		self.surf_mean = mean(self.surf,axis=0)

#---Undulation spectra functions

	def calculate_undulations(self,removeavg=0,redundant=1,whichframes=None,qmagfilter=None,
		fitbest=False,fitlims=None,forcekappa=True,peri=False):
		'''Fourier transform surface heights.'''
		if fitlims == None: 
			fitbest = False
			fitlims = [1024,8]
		if whichframes == None:
			framerange = range(len(self.surf))
		else:
			framerange = whichframes
		lenscale = self.lenscale
		status('status: lenscale = '+str(lenscale))
		#---undulations
		self.undulate_raw = []
		if removeavg == 0:
			for k in framerange:
				status('status: Fourier transform, frame ',i=k,looplen=len(framerange))
				if redundant == 1:
					self.undulate_raw.append(fft.fftshift(fft.fft2(array(self.surf[k])[:-1,:-1]/lenscale)))
				elif redundant == 2:
					self.undulate_raw.append(fft.fftshift(fft.fft2(array(self.surf[k])[:-2,:-2]/lenscale)))
				else:
					self.undulate_raw.append(fft.fftshift(fft.fft2(array(self.surf[k])/lenscale)))
		else:
			self.calculate_average_surface()
			for k in range(len(self.surf)):
				status('status: Fourier transform, frame ',i=k,looplen=len(framerange))
				if redundant == 1:
					relative = array(self.surf[k])/lenscale-array(self.surf_mean)/lenscale
					self.undulate_raw.append(fft.fftshift(fft.fft2(relative[:-1,:-1])))
				elif redundant == 2:
					relative = array(self.surf[k])/lenscale-array(self.surf_mean)/lenscale
					self.undulate_raw.append(fft.fftshift(fft.fft2(relative[:-2,:-2])))
				else:
					self.undulate_raw.append(fft.fftshift(fft.fft2(array(self.surf[k])/lenscale-
						array(self.surf_mean)/lenscale)))
		#---peristalsis
		if peri:
			self.undulate_peri_raw = []
			for k in framerange:
				status('status: Fourier transform, frame ',i=k,looplen=len(framerange))
				if redundant == 1:
					self.undulate_peri_raw.append(fft.fftshift(fft.fft2(
						array(self.surf_thick[k])[:-1,:-1]/lenscale)))
				elif redundant == 2:
					self.undulate_peri_raw.append(fft.fftshift(fft.fft2(array(
						self.surf_thick[k])[:-2,:-2]/lenscale)))
				else:
					self.undulate_peri_raw.append(fft.fftshift(fft.fft2(array(
						self.surf_thick[k])/lenscale)))
		#---calculate 2D spectrum
		m,n = shape(self.undulate_raw)[1:]
		#---follows scipy DFFT convention on even/odd location of Nyquist component
		spec_center = [int(i/2) for i in shape(self.undulate_raw)[1:]]
		lxy = array([self.vec(i) for i in self.surf_index])/lenscale
		yv,xv = meshgrid(range(n),range(m))
		qs = [(sqrt((2*pi*(array(xv)-spec_center[0])/lxy[f][0])**2+
			(2*pi*(array(yv)-spec_center[1])/lxy[f][1])**2)) for f in range(len(lxy))]
		self.undulate_hqhq2d = mean(array(1.*(abs(array(self.undulate_raw))/double(m*n))**2),axis=0)
		self.undulate_qmag2d = mean(qs,axis=0)
		self.undulate_spec1d = array([self.undulate_qmag2d,self.undulate_hqhq2d]).T.reshape(m*n,2)
		if peri:
			self.undulate_peri_hqhq2d = mean(array(1.*(abs(array(self.undulate_peri_raw))/double(m*n))**2),
				axis=0)
			self.undulate_peri_qmag2d = mean(qs,axis=0)
			self.undulate_peri_spec1d = array([self.undulate_peri_qmag2d,
				self.undulate_peri_hqhq2d]).T.reshape(m*n,2)
		spec1d = array([i for i in array(self.undulate_spec1d) if i[0] != 0.])
		specsort = spec1d[np.lexsort((spec1d[:,1],spec1d[:,0]))]
		#---fitting the best points based on prescribed limits
		if fitbest and qmagfilter == None:	
			best_rmsd,best_endpost = 10**10,0
			for endpost in range(0,len(specsort)/fitlims[0]):
				specfilter = specsort[endpost:len(specsort)/fitlims[1]]
				[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
				if sum([(exp(az)*(specfilter[i,0]**bz)-specfilter[i,1])**2 
					for i in range(len(specfilter))])/len(specfilter[:,0]) < best_rmsd:
					best_rmsd = sum([(exp(az)*(specfilter[i,0]**bz)-specfilter[i,1])**2 
						for i in range(len(specfilter))])/len(specfilter[:,0])
					best_endpost = endpost
			qmagfilter = [specsort[best_endpost,0],specsort[len(specsort)/fitlims[1],0]]
		#---otherwise set the qmagfilter based on the number of points according to fitlims
		elif not fitbest and qmagfilter == None:
			qmagfilter = [10**-10,specsort[len(specsort)/fitlims[1],0]]
		specfilter = array(filter(lambda x: x[0] >= qmagfilter[0] 
			and x[0] <= qmagfilter[1],self.undulate_spec1d))
		area = double(mean([self.vec(i)[0]*self.vec(i)[1] for i in self.surf_index])/lenscale**2)
		[bz,az]=numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
		leftcom = [mean(log(specfilter[:,0])),mean(log(specfilter[:,1]))]
		az_enforced = leftcom[1]+4.*leftcom[0]
		self.undulate_qmagfilter = qmagfilter
		if forcekappa:
			self.undulate_kappa = 1/exp(az_enforced)/area
		else:
			self.undulate_kappa = 1/exp(az)/area

#---Lipid packing and tilt properties

	def tilter(self,selector,director,framecount=None,start=None,end=None,
		skip=None,protein_selection=None,residues=None):
		'''UNDER CONSTRUCTION: Lipid tilt calculator.'''
		if framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		#self.griddims = [int(round(self.vec(1)[0]/rounder)),int(round(self.vec(1)[1]/rounder))]
		result_data = MembraneData('tilt')
		result_data_position = MembraneData('lipid_positions')
		for fr in range(start,end,skip):
			print 'status: calculating surface normals and tilt for frame: '+str(fr)
			starttime = time.time()
			self.gotoframe(fr)
			#self.calculate_midplane(selector,k,rounder=rounder,interp=interp,residues=residues)
			#---start mod
			starttime = time.time()
			vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
			find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices 
				if x in simplex for indx in simplex if indx !=x))
			print 'status: getting points'
			topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for i in self.monolayer_residues[0]])
			#botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
			#	for i in self.monolayer_residues[1]])
			toptailxyz = array([mean(self.universe.residues[i].selectAtoms(''.join([i+' or ' 
				for i in director[1:-1]]+[director[-1]])).coordinates(),axis=0) 
				for i in self.monolayer_residues[0]])
			topxyz_wrapped = self.wrappbc(topxyz,self.vec(fr),mode='grow')
			dtri = scipy.spatial.Delaunay(topxyz_wrapped[:,0:2])
			pts = topxyz_wrapped
			point_permute = [[(i+j)%3 for i in range(3)] for j in range(3)]
			print 'status: calculating triangle face normals'
			trifaces = [[np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]) 
				for i in point_permute] for j in dtri.simplices]
			print 'status: calculating simplex areas'
			simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
			print 1./60.*(time.time()-starttime)	
			ptsareas = np.zeros(len(dtri.points))
			ptsnorms = np.zeros([len(dtri.points),3])
			print 'status: summing'
			for s in range(len(dtri.simplices)):
				for p in range(3):
					ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
						if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
				ptsareas[dtri.simplices[s][p]] += simp_areas[s]
			ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
			print 'status: calculating angle'
			vecslipids = [toptailxyz[i]-topxyz[i] for i in range(len(topxyz))]
			#print 1./60.*(time.time()-starttime)
			#meshpoints(array(ptsnorms)+array(pts))
			#meshplot(topxyz_wrapped)
			angles = [1./pi*arccos(np.dot(vecslipids[i],ptsnorms[i])/np.linalg.norm(
				vecslipids[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
			#plotthat = [[topxyz[i][0],topxyz[i][1],50*angles[i]] for i in range(len(topxyz))]
			#meshplot(plotthat)
			#areas = [1./3.*sum() for i in range(len(topxyz))]
			[1./pi*arccos(np.dot(vecslipids[i],ptsnorms[i])/np.linalg.norm(
				vecslipids[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
			result_data.add([[[angles],[]],[[ptsareas],[]]],[fr])
			result_data_position.add([topxyz,[]],[fr])
			if protein_selection != None:
				self.protein.append(self.universe.selectAtoms(protein_selection).coordinates())
				self.protein_index.append(fr)
		self.store.append(result_data)
		self.store.append(result_data_position)
		del result_data
		del result_data_position
		#---end mod

#---Radial distribution functions
	
	def batch_gr_lipid_ion(self,selector,start=None,end=None,skip=None,framecount=None,label='',mode=None,
		monolayer_rep=None,monos=None):
		'''POSSIBLY DEPRECATED: Calculate the radial distribution function between lipids and ions.'''
		if framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		print 'status: starting g(r) calculation'
		if end == 0 or end == None: end = self.nframes
		#---Identify monolayers
		if monolayer_rep != None: self.monolayer_rep = monolayer_rep
		if self.monolayer_residues == []:
			self.identify_monolayers()
		#---Select atoms
		allselect_lipids = self.universe.selectAtoms(selector[0])
		allselect_ions = self.universe.selectAtoms(selector[1])
		validresids = list(set.intersection(set(self.monolayer_residues[0]),
			set([i-1 for i in allselect_lipids.resids()])))
		sel1 = sum([allselect_lipids.residues[list(allselect_lipids.resids()).index(i+1)].\
			selectAtoms(selector[0]) for i in validresids])
		validresids = list(set.intersection(set(self.monolayer_residues[1]),
			set([i-1 for i in allselect_lipids.resids()])))
		sel2 = sum([allselect_lipids.residues[list(allselect_lipids.resids()).index(i+1)].\
			selectAtoms(selector[0]) for i in validresids])
		self.selections.append(sel1)
		self.selections.append(sel2)
		self.selections.append(allselect_ions)
		#---Process the frames
		result_data = MembraneData('grvoronoi' if mode == 'voronoi_bin' else 'gr',label=label)
		for k in range(start,end,skip):
			print 'status: calculating RDF, frame: '+str(k)
			if monos == None:
				result_data.add([self.calculate_gr_lipid_ion_1d(k,whichmono=0,detectside=1,mode=mode),
					self.calculate_gr_lipid_ion_1d(k,whichmono=1,detectside=1,mode=mode)],[k])
			else:
				result_data.add([self.calculate_gr_lipid_ion_1d(k,whichmono=mono,detectside=1,mode=mode) 
					for mono in [0,1] if mono in monos],[k])
		self.store.append(result_data)
		del result_data
		#---Clear selections
		self.selections = []
	
	def calculate_gr_lipid_ion_1d(self,frameno,whichmono=-1,detectside=0,
		pbcadjust=1,mode=None,duplicate=True):
		'''POSSIBLY DEPRECATED: Framewise radial distribution calculator.'''
		if (whichmono == 0 and self.selections[0] == 0) or (whichmono == 1 and self.selections[1] == 0):
			return [[],[]]
		if detectside != 0:
			if whichmono == 0:
				lipids_in = array(self.get_points(frameno,selection_index=0))
				ions_in_all = self.get_points(frameno,selection_index=2)
				bilayer_z = self.locate_bilayer(frameno)
				if duplicate == True:
					ions_in = array([c+[0,0,1*(((c[2]<bilayer_z)*self.vec(frameno)[2])-0*bilayer_z)] 
						for c in ions_in_all])
				else:
					ions_in = [c+[0,0,-1*(((int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*\
						self.vec(frameno)[2])-0*bilayer_z)] for c in ions_in_all if 
						((c[2]+-1*(int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*\
						self.vec(frameno)[2])-bilayer_z) > 0.0]
			elif whichmono == 1:
				lipids_in = array(self.get_points(frameno,selection_index=1))
				ions_in_all = self.get_points(frameno,selection_index=2)
				bilayer_z = self.locate_bilayer(frameno)
				if duplicate == True:
					ions_in = array([c+[0,0,-1*(((c[2]>bilayer_z)*self.vec(frameno)[2])-0*bilayer_z)] 
						for c in ions_in_all])
				else:				
					ions_in = [c+[0,0,-1*(((int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*\
						self.vec(frameno)[2])-0*bilayer_z)] for c in ions_in_all if \
						((c[2]+-1*(int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*\
						self.vec(frameno)[2])-bilayer_z) <= 0.0]
		if mode == None:
			subject = array(ions_in)[:,2:3]
			snap = subject[:]
			pairdistframe = []
			distancesout = scipy.spatial.distance.cdist(lipids_in[:,2:3],snap)
			return distancesout
		elif mode == 'voronoi_bin':
			tree = scipy.spatial.KDTree(lipids_in[:,0:2])
			tmp = [tree.query(i)[1] for i in array(ions_in)[:,0:2]]
			vordists = [tree.query(i)[0] for i in array(ions_in)[:,0:2]]
			binned = [[i for i,x in enumerate(tmp) if x == b] for b in range(len(lipids_in))]
			pdistsz = [[abs(lipids_in[i][2]-ions_in[binned[i][j]][2]) for j in range(len(binned[i]))]
				for i in range(len(lipids_in)) if len(binned[i]) != 0]
			pdists2d = [[linalg.norm(lipids_in[i][0:2]-ions_in[binned[i][j]][0:2]) 
				for j in range(len(binned[i]))]
				for i in range(len(lipids_in))  if len(binned[i]) != 0]
			return [pdistsz,pdists2d]
			
	def batch_gr_lipid(self,selector,start=None,end=None,skip=None,framecount=None,label='',mode=None,
		monolayer_rep=None,monos=None):
		'''POSSIBLY DEPRECATED: Calculate the radial distribution function between lipids in 2D.'''
		if framecount == None:
			if end == None: end = self.nframes
			if start == None: start = 0
			if skip == None: skip = 1
		else:
			start = 0
			end = self.nframes
			skip = int(float(self.nframes)/framecount)
			skip = 1 if skip < 1 else skip
		print 'status: starting g(r) calculation'
		if end == 0 or end == None: end = self.nframes
		#---Identify monolayers
		if monolayer_rep != None: self.monolayer_rep = monolayer_rep
		if self.monolayer_residues == []:
			self.identify_monolayers()
		#---Select atoms
		allselect_lipids1 = self.universe.selectAtoms(selector[0])
		validresids = list(set.intersection(set(self.monolayer_residues[0]),
			set([i-1 for i in allselect_lipids1.resids()])))
		sel1 = sum([allselect_lipids1.residues[allselect_lipids1.resids().index(i+1)].selectAtoms(selector[0])
			for i in validresids])
		validresids = list(set.intersection(set(self.monolayer_residues[1]),
			set([i-1 for i in allselect_lipids1.resids()])))
		sel2 = sum([allselect_lipids.residues[allselect_lipids.resids().index(i+1)].selectAtoms(selector[0])
			for i in validresids])
		self.selections.append(sel1)
		self.selections.append(sel2)
		#---finish code in jot-lipid-gr.py and integrate it here
		allselect_lipids2 = self.universe.selectAtoms(selector[1])
		self.selections.append(allselect_ions)
		#---Process the frames
		result_data = MembraneData('grvoronoi' if mode == 'voronoi_bin' else 'gr',label=label)
		for k in range(start,end,skip):
			print 'status: calculating RDF, frame: '+str(k)
			if monos == None:
				result_data.add([self.calculate_gr_lipid_ion_1d(k,whichmono=0,detectside=1,mode=mode),
					self.calculate_gr_lipid_ion_1d(k,whichmono=1,detectside=1,mode=mode)],[k])
			else:
				result_data.add([self.calculate_gr_lipid_ion_1d(k,whichmono=mono,detectside=1,mode=mode) 
					for mono in [0,1] if mono in monos],[k])
		self.store.append(result_data)
		del result_data
		#---Clear selections
		self.selections = []

#---Making regular, triangulated meshes

	def wrappbc(self,points,vecs,dims=[0,1],mode=None,growsize=0.2):
		'''Adjusts input points to wrap or reflect over periodic boundaries in multiple ways.'''
		#---Takes anything outside the box, and adds it to the other side. Useful for non-rectangular meshes.
		if mode == 'add_oversteps':
			ans = []
			for p in points:
				if sum([(int(p[i]/vecs[i])%2) for i in dims]) > 0:
					ans.append([p[i]-((int(p[i]/vecs[i])%2)*vecs[i] if i in dims else 0) \
						for i in range(shape(points)[1])])
			return concatenate((points,array(ans)))
		#---Increases system size by the specified growsize percentage.
		elif mode == 'grow':
			ans = []
			over = [[0,1],[1,0],[-1,0],[0,-1],[-1,1],[1,-1],[1,1],[-1,-1]]
			for p in points:
				for move in over:
					ptmove = list([p[i]+vecs[i]*move[i] if i in dims else p[i] for i in range(3)])
					if ((ptmove[0] >= -vecs[0]*growsize) and (ptmove[0] <= vecs[0]*(1.+growsize)) \
						and (ptmove[1] >= -vecs[1]*growsize) and (ptmove[1] <= vecs[1]*(1.+growsize))):
						ans.append(ptmove)
			return concatenate((points,array(ans)))
		#---Includes adjacent (but not corner) meshes.
		elif mode == 'embiggen':
			ans = []
			for p in points:
				for tr in [[1,0],[0,1],[-1,0],[0,-1]]:
					ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(3)])
			return concatenate((points,array(ans)))
		elif mode == 'nine':
			ans = []
			for p in points:
				for tr in [[1,0],[0,1],[-1,0],[0,-1],[-1,1],[-1,-1],[1,1],[1,-1]]:
					ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(3)])
			return concatenate((points,array(ans)))
		else:
			return array([\
				[p[i]-((int(p[i]/vecs[i])%2)*vecs[i] if i in dims else 0) \
				for i in range(shape(points)[1])] for p in points])
				
	def binner(self,points,vecs,grid):
		'''Binning function to reduce the complexity of high-resolution 3-space points.'''
		ti1 = linspace(0,vecs[0],grid[0])
		ti2 = linspace(0,vecs[1],grid[1])
		rounder_vecs = [ti1[1],ti2[1]]
		binner = [[[] for j in range(grid[1])] for i in range(grid[0])]
		for pt in points:
			binner[int(round(pt[0]/rounder_vecs[0]))][int(round(pt[1]/rounder_vecs[1]))].append(pt[2])
		binnedpts = []
		for i in range(grid[0]):
			for j in range(grid[1]):
				if binner[i][j] != []:
					binnedpts.append([ti1[i],ti2[j],mean(binner[i][j])])
		return array(binnedpts)
		
	def makemesh(self,data,vecs,grid,method='best',fine=None):
		'''Approximates an unstructured mesh with a regular one.'''
		if method == 'best': method = 'bilinear_cubic'
		if method == 'bilinear_triangle_interpolation':
			dat1 = self.wrappbc(data,vecs,mode='grow')
			dat2 = Delaunay(dat1[:,0:2])
			xypts = array([[i,j] for i in linspace(0,vecs[0],grid[0]) for j in linspace(0,vecs[1],grid[1])])
			results = []
			for i in range(len(xypts)):
				print 'status: interpolating point: '+str(i)
				dat2 = Delaunay(dat1[:,0:2])
				mytripts = dat1[dat2.simplices[dat2.find_simplex(xypts[i])]]
				height = scipy.interpolate.griddata(mytripts[:,0:2],mytripts[:,2],[xypts[i]],method='linear')
				results.append([xypts[i][0],xypts[i][1],height])
			return array(results)
		elif method == 'bilinear_triangle_interpolation_manual':
			#---Manual bilinear interpolation. Same speed as griddata in method 1 (above).
			dat1 = self.wrappbc(data,vecs,mode='grow')
			dat2 = Delaunay(dat1[:,0:2])
			xypts = array([[i,j] for i in linspace(0,vecs[0],grid[0]) for j in linspace(0,vecs[1],grid[1])])
			results = []
			for i in range(len(xypts)):
				print 'status: interpolating point: '+str(i)
				dat2 = Delaunay(dat1[:,0:2])
				t = dat1[dat2.simplices[dat2.find_simplex(xypts[i])]]
				det = t[0][0]*t[1][1]-t[1][0]*t[0][1]+t[1][0]*t[2][1]-t[2][0]*t[1][1]+t[2][0]*t[0][1]\
					-t[0][0]*t[2][1]
				a = (((t[1][1]-t[2][1])*t[0][2]+(t[2][1]-t[0][1])*t[1][2]+(t[0][1]-t[1][1])*t[2][2]))/det
				b = (((t[2][0]-t[1][0])*t[0][2]+(t[0][0]-t[2][0])*t[1][2]+(t[1][0]-t[0][0])*t[2][2]))/det 
				c = (((t[1][0]*t[2][1]-t[2][0]*t[1][1])*t[0][2]+(t[2][0]*t[0][1]-t[0][0]*t[2][1])*t[1][2]\
					+(t[0][0]*t[1][1]-t[1][0]*t[0][1])*t[2][2]))/det
				height = a*xypts[i][0]+b*xypts[i][1]+c
				results.append([xypts[i][0],xypts[i][1],height])
			return array(results)
		elif method == 'bilinear':
			starttime = time.time()
			xypts = array([[i,j] for i in linspace(0,vecs[0],grid[0]) for j in linspace(0,vecs[1],grid[1])])
			interp = scipy.interpolate.LinearNDInterpolator(data[:,0:2],data[:,2],fill_value=0.0)
			return array([[i[0],i[1],interp(i[0],i[1])] for i in xypts])
		elif method == 'bilinear_cubic':
			starttime = time.time()
			xypts = array([[i,j] for i in linspace(0,vecs[0],grid[0]) for j in linspace(0,vecs[1],grid[1])])
			interp = scipy.interpolate.LinearNDInterpolator(data[:,0:2],data[:,2],fill_value=0.0)
			bilinear_pts = array([[i[0],i[1],interp(i[0],i[1])] for i in xypts])
			result = scipy.interpolate.griddata(bilinear_pts[:,0:2],bilinear_pts[:,2],bilinear_pts[:,0:2],
				method='cubic')
			return array([[bilinear_pts[i,0],bilinear_pts[i,1],result[i]] for i in range(len(result))])
		elif method == 'rbf':
			subj = self.wrappbc(data,vecs,mode='grow')
			print 'status: generating radial basis function object'
			rbfobj = scipy.interpolate.Rbf(subj[:,0],subj[:,1],subj[:,2],epsilon=1.2,function='gaussian')
			ti1 = linspace(0,vecs[0],grid[0]);ti2 = linspace(0,vecs[1],grid[1])
			XI, YI = meshgrid(ti1, ti2)
			print 'status: surfacing on a regular grid'
			ZI = rbfobj(XI,YI)
			return self.unzipgrid(transpose(ZI),vecs=vecs,reverse=0)
		elif method == 'sharp_cubic_spline':
			'''
			Added this method on 2014.08.14 in order to dramatically increase the resolution of the mesh.
			Currently the sharpness parameter is hard-coded but this should eventually become an argument.
			'''
			if fine == None: fineness = 500*1j,500*1j
			else: fine = fine*1j,fine*1j
			
			gridx,gridy = mgrid[0:vecs[0]:fine[0],0:vecs[1]:fine[1]]
			points = data[:,:2]
			values = data[:,2]
			heights = scipy.interpolate.griddata(points,values,(gridx,gridy),method='cubic')
			unzip = self.unzipgrid(transpose(heights),vecs=vecs,reverse=0)
			return unzip
		
#---Transforming points in 3-space to heights in a grid, and vis-versa

	def unzipgrid(self,surf,vecs=None,grid=None,rounder_vecs=[1.0,1.0],reverse=0):
		'''Turns a 2D array into a set of points in 3-space.'''
		if type(surf) != ndarray:
			surf = array(surf)
		grid = [shape(surf)[i] for i in range(2)]
		if reverse != 0: grid = grid[::-1];
		if vecs != None:
			rounder_vecs = [vecs[i]/(grid[i]-1) for i in range(2)]
		replotin = surf
		surfreplot = []
		for i in range(grid[0]):
				for j in range(grid[1]):
					surfreplot.append([i*rounder_vecs[0],j*rounder_vecs[1],replotin[i,j]])
		surfreplot = array(surfreplot)
		return surfreplot
		
	def rezipgrid(self,xyz,vecs=None,frameno=0,grid=None,rounder_vecs=[1.0,1.0],
		reverse=0,diff=False,whichind=2):
		'''Turns a regular set of points in 3-space into a 2D matrix.'''
		#---Nb this is the source of the transpose error, which needs fixed.
		#---Nb it looks like I fixed this but didn't write it down until afterwards.
		#---Modifications in the following section for compatibility with the general interpolation function.
		#---Nb "diff" describes whether we are handling redundant points. Needs explained.
		#---Nb I removed all references to diff in the process of fixing script-meso-coupling code.
		if grid == None and diff != True: grid = self.griddims
		elif grid == None and diff == True: grid = self.griddims[0],self.griddims[1]
		if vecs == None: vecs = self.vec(frameno)
		steps = ([vecs[i]/(grid[i]-1) for i in range(2)] 
			if diff == False else [vecs[i]/(grid[i]) for i in range(2)])
		poslookup = [[xyz[i][j]/steps[j] for j in range(2)] for i in range(len(xyz))]
		surfgrid = [[0. for i in range(grid[1])] for j in range(grid[0])]
		for i in range(len(xyz)): 
			#---Note: added the round command below changes calculate_midplane time from 0.25 to 0.35!
			if int(poslookup[i][0]) < grid[0] and int(poslookup[i][1]) < grid[1]:
				surfgrid[int(round(poslookup[i][0]))][int(round(poslookup[i][1]))] = xyz[i][whichind]
		return surfgrid
		
#---Retrieving data

	def getdata(self,label):
		return self.store[[i.description for i in self.store].index(label)]
		
	def avail(self):
		return [i.description for i in self.store]

#---EXTERNAL FUNCTIONS, PICKLE THE DATA
#-------------------------------------------------------------------------------------------------------------

def pickledump(obj,filename,directory=''):
	'''Pickles an object to a text file.'''
	if obj.__class__.__name__ == 'MembraneSet':
		obj.picklename = filename.strip('pkl.').strip('.pkl')
		obj.universe = []
		obj.xyzs = []
	if os.path.isfile(directory+filename):
		for i in range(1,100):
			latestfile = directory+'#'+filename+'.'+('%02d' % i)+'#'
			if not os.path.isfile(directory+'#'+filename+'.'+('%02d' % i)+'#'):
				os.rename(directory+filename,latestfile)
				print 'status: backing up a pre-existing pickle file with the same name to number '+str(i)
				break
	#---clear bulky data automatically before writing
	if hasattr(obj,'universe') == True: obj.universe = []
	if hasattr(obj,'xyzs') == True: obj.xyzs = []
	fp = open(directory+filename, 'w')
	pickle.dump(obj,fp)
	fp.close()

def unpickle(filename):
	'''Un-pickles an object from a text file.'''
	status('status: seeking or loading '+filename)
	if os.path.isfile(filename):
		fp = open(filename, 'r')
		x = pickle.load(fp)
		if hasattr(x,'picklename') == False and \
			type(x) != list and type(x) != ndarray and type(x) != dict:
			x.picklename = filename.strip('pkl.').strip('.pkl')
		elif type(x) != list and type(x) != ndarray and type(x) != dict:
			if x.picklename == '': 
				x.picklename = filename.strip('pkl.').strip('.pkl')
		fp.close()
		return x
	else:
		status('status: file does not exist')
		return
		
def binarysave(obj,filename,directory=''):
	'''Pickles an object to a text file.'''
	if os.path.isfile(directory+filename):
		for i in range(1,100):
			latestfile = directory+'#'+filename+'.'+('%02d' % i)+'#'
			if not os.path.isfile(directory+'#'+filename+'.'+('%02d' % i)+'#'):
				os.rename(directory+filename,latestfile)
				print 'status: backing up a pre-existing h5 binary file with the same name to number '+str(i)
				break
	h5f = h5py.File(directory+filename,'w')
	h5f.create_dataset('data',data=obj)
	h5f.close()

def unbinary(filename):
	'''Un-pickles an object from a text file.'''
	status('status: seeking or loading '+filename)
	if os.path.isfile(filename):
		h5f = h5py.File(filename,'r')
		data = h5f['data'][:]
		h5f.close()
		return data
	else:
		status('status: file does not exist')
		return
	

#!/usr/bin/env python

import numpy
import collections
import itertools

#---MEMBRANEDATA CLASS
#-------------------------------------------------------------------------------------------------------------

def flatten(x,loops=None):
	'''Flatten an arbitrarily deep list. Necessary for arbitrary slicing in MembraneData.'''
	if isinstance(x, collections.Iterable):
		return [a for i in x for a in flatten(i)]
	else:
		return [x]

class MembraneData:
	'''A class which holds a membrane simulation.'''
	def __init__(self,name,label=None):
		#---Descriptors
		self.calctype = name
		self.pickle_name = ''
		self.system_name = ''
		if self.calctype == 'grvoronoi':
			self.struct = {'frame':0,'monolayer':1,'direction':2,'lipid':3}
			self.struct_opts = {'direction' : {'z':0,'2d':1}}
		elif self.calctype == 'gr':
			self.struct = {'frame':0,'monolayer':1}
		elif self.calctype == 'triangles':
			self.struct = {'frame':0,'monolayer':1,'type':2}
			self.struct_opts = {'type' : {'points':0,'lines':1,'n_perfect_triangles':2}}
		elif self.calctype == 'cells':
			self.struct = {'frame':0,'monolayer':1,'type':2}
			self.struct_opts = {'type' : {'points':0,'voronoi':1,'areas':2}}
		#---retire the following data type, representing a verbatim port from old method
		elif self.calctype == 'dimple_filter':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'points':0,'domain':1,'protptsplus':2,'domainw':3,'frameno':4}}
		#---retired method, formerly named "dimple" so be careful not to confuse them
		elif self.calctype == 'dimple_for':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}}
		elif self.calctype == 'tilt_deprecated':
			self.struct = {'frame':0,'type':1,'monolayer':2,'lipid':3}
			self.struct_opts = {'type' : {'angle':0,'area':1}}
		elif self.calctype == 'tilefilter_area_v1':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'positive':0,'negative':1}}
		#---stores surface normal vectors
		elif self.calctype == 'surfnorms':
			self.struct = {'frame':0,'monolayer':1,'type':2,'lipid':3}
			self.struct_opts = {'type':{'normals','areas'}}
		elif self.calctype == 'tilts':
			self.struct = {'frame':0,'tail':1,'monolayer':2,'lipid':3}
			self.struct_opts = {'tail' : {'a':0,'b':1}}
		elif self.calctype == 'protdists':
			self.struct = {'frame':0,'monolayer':2,'lipid':3}
		elif self.calctype == 'lipid_positions':
			self.struct = {'frame':0,'type':1,'monolayer':2,'lipid':3}
			self.struct_opts = {'type' : {'position':0,'mindist':1}}
		elif self.calctype == 'dimple':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}}
		elif self.calctype == 'gr2d':
			self.struct = {'frame':0,'monolayer':1,'lipid':2}
		elif self.calctype == 'xyz_assoc':
			self.struct = {'frame':0,'monolayer':1,'lipid':2}
		#---type: c0map holds data in stressmap pickles after integrating the voxel-wise stress tensors
		elif self.calctype == 'c0map':
			self.struct = {'frame':0}
		#---type: topography_transform holds an array of minimum distance vs surface point height
		#---Nb this is an alternative to the topographycorrelate class
		elif self.calctype == 'topography_transform':
			self.struct = {'frame':0}
		#---type: ionskate for storing ice-skating ion MSD decompositions
		#---Nb the mastermsd_zones array is very jagged so easy lookups might not work
		#---Nb this was basically scrapped due to memory issues
		elif self.calctype == 'ionskate':
			self.struct = {'type':0,'zone':1,'ion':2,'deltat':3,'start_frame':4}
			self.struct_opts = {'type': {'mastermsd_zones':0,'distsxy':1,'distsz':2}}
		elif self.calctype == 'spanangle':
			self.struct = {'frame':0,'resid':1,'headspan':2,'headangle':3}
		elif self.calctype == 'spanangle2':
			self.struct = {'frame':0,'resid':1,'type':2}
			self.struct_opts = {'type': {'headspan':0,'headangle':1,'time':2,'resid':3}}
		#---type: updated/advanced dimple fitting (dimple2 pkl objects)
		elif self.calctype == 'dimple2':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}}
		#---type: updated/advanced dimple fitting (dimple3 pkl objects)
		elif self.calctype == 'dimple3':
			self.struct = {'frame':0,'type':1}
			self.struct_opts = {'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}}
		#---type: updated object which holds height distributions near certain neighborhoods (topography)
		elif self.calctype == 'topography3':
			self.struct = {'frame':0,'heights':1}
		#---default values for label if no extra descriptors
		if label != None:
			self.description = label 
		else:
			self.description = name
		#---Data
		self.data = []
		self.label = []
		self.notes = []
	#---Add a piece of data
	def add(self,item,descriptor):
		self.data.append(item)
		self.label.append(descriptor)
	#---Add a note
	def addnote(self,notable):
		self.notes.append(notable)
	#---retrieve a note
	def getnote(self,note):
		availnotes = [i[0] for i in self.notes]
		if note in availnotes:
			return self.notes[availnotes.index(note)][1]
		#else: print 'missing note'
	#---Retrieve data
	def get(self,*args,**kwargs):
		#---Flatten the data as necessary, or depending on the calculation
		if 'flat' in kwargs.keys():
			flat = kwargs['flat']
		else:
			if (self.calctype == 'grvoronoi' or self.calctype == 'gr'):
				flat = True
			elif (self.calctype == 'triangles' or self.calctype == 'cells'):
				flat = False
			else:
				flat = False
		#---Identify dimensions to slice, and what their slice will be
		keydims = []
		slices = []
		structs = []
		#---Handle fluid arguments
		if numpy.shape(args) == (2,2):
			pass
		elif len(numpy.shape(args)) == 3 and numpy.shape(args)[0] == 1:
			args = tuple(args[0])
		elif numpy.shape(args)[1] > 2:
			args = tuple([[args[0][j],args[0][j+1]] for j in range(0,len(args[0]),2)])
		for spec in args:
			keydims.append(self.struct[spec[0]])
			structs.append(spec[0])
			slices.append(spec[1])
		#---Generate a slice variable
		slicer = []
		for i in range(len(numpy.shape(self.data))):
			if i in keydims:
				if type(slices[keydims.index(i)]) != str:
					slicer.append(slices[keydims.index(i)])
				else:
					slicer.append(self.struct_opts[structs[keydims.index(i)]][slices[keydims.index(i)]])
			else:
				slicer.append(slice(0,None))
		slicer = tuple(slicer)
		#---Return the data, flattened if desired
		if flat:
			return flatten(numpy.array(self.data)[slicer])
		else:
			return numpy.array(self.data)[slicer]


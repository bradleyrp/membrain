#!/usr/bin/env python

import numpy
import collections
import itertools

#---master dictionary which provides specifications for different data types
#---note that the struct and struct_opts let you specify various dimensions for your numpy array result
#---note that this slicing makes it possible to reference data by the name of a particular dimension
#---? include more helpful explanation of how slicing works
master_datatypes = {
	'gr2d':{
		'struct':{'frame':0,'monolayer':1,'lipid':2},
		'struct_opts':None,
		'description':
		'Basic object which actually contains 3D radial distribution function (i.e. g(r)) data for '+\
		'AAMD bilayers. The root data is a normalized histogram or probability curve of the likelihood of '+\
		'reaching another particle. Particle specs are in the pkl name and the notes.',
		},
	'cells':{
		'struct':{'frame':0,'monolayer':1,'type':2},
		'struct_opts':{'type':{'points':0,'voronoi':1,'areas':2}},
		'description':'Holds the list of raw points, the Voronoi mesh, and the area per neighborhood. '+\
		'These are often stored in separate pickles, with the other fields blanked out.',
		},
	'spanangle2':{
		'struct':{'frame':0,'resid':1,'type':2},
		'struct_opts':{'type': {'headspan':0,'headangle':1,'time':2,'resid':3}},
		'description':'Holds the headgroup angle and headgroup span for use in AAMD bilayer simulations',
		},
	'c0map':{
		'struct':{'frame':0},
		'struct_opts':None,
		'description':'holds interpolated C_0 map from mesoscale simulations',
		},
	'surfnorms':{
		'struct':{'frame':0,'monolayer':1,'type':2,'lipid':3},
		'struct_opts':{'type':{'normals','areas'}},
		'description':'stores surface normal vectors',
		},
	'collect_c0maps':{
		'struct':{'frame':0},
		'struct_opts':None,
		'description':'c0map holds data in stressmap pickles after integrating the voxel-wise stress tensors',
		},
	'polyangles':{
		'struct':{'frame':0,'angles':1},
		'struct_opts':None,
		'description':'stores the angles along a polymer chain',
		},
    'bridge':{
		'struct':{'frame':0,'monolayer':1,'lipid':2},
		'struct_opts':None,
		'description':'Holds the output of ion-lipid bridge calculations. Particle specs are in the pkl '+\
		'name and the notes.',
		},
	}

#---MEMBRANEDATA CLASS
#-------------------------------------------------------------------------------------------------------------

def flatten(x,loops=None):
	'''Flatten an arbitrarily deep list. Necessary for arbitrary slicing in MembraneData.'''
	if isinstance(x,collections.Iterable):
		return [a for i in x for a in flatten(i)]
	else:
		return [x]

class MembraneData:
	'''A class which holds a membrane simulation or associated data.'''
	def __init__(self,name,label=None):
		#---standard descriptors
		self.calctype = name
		self.pickle_name = ''
		self.system_name = ''
		#---look up the type definitions from the master list
		if name in master_datatypes.keys():
			self.struct = (master_datatypes[name])['struct']
			self.struct_opts = (master_datatypes[name])['struct_opts']
		else: raise Exception('except: data type does not exist')
		#---an extra descriptor for the object, distinct from label, which labels added data pieces
		if label != None: self.description = label 
		else: self.description = name
		#---data objects
		self.data = []
		self.label = []
		self.notes = []
	#---add a piece of data with a label, usually the frame number as an extra check
	def add(self,item,descriptor):
		self.data.append(item)
		self.label.append(descriptor)
	#---add a note, usually information from analysis_descriptors about the corresponding simulation
	def addnote(self,notable):
		self.notes.append(notable)
	#---retrieve a note
	def getnote(self,note):
		availnotes = [i[0] for i in self.notes]
		if note in availnotes: return self.notes[availnotes.index(note)][1]
		else: return None
	#---retrieve data
	#---note that this function encompasses the entire utility of membraindata
	#---...in which the membraindata object provides a secure way to recall dimensions of the incoming 
	#---...data so that the user never gets confused about the physical meaning of each dimension
	def get(self,*args,**kwargs):
		#---flatten some types of outgoing data
		#---? move this parameter to the dictionary
		if 'flat' in kwargs.keys():	flat = kwargs['flat']
		else:
			if (self.calctype == 'grvoronoi' or self.calctype == 'gr'): flat = True
			elif (self.calctype == 'triangles' or self.calctype == 'cells'): flat = False
			else: flat = False
		#---identify dimensions to slice, and what their slice will be
		#---slice referes to the actual slice along the dimension specified by struct in the dictionary
		keydims = []
		slices = []
		structs = []
		#---handle a few different types of arguments
		if numpy.shape(args) == (2,2): pass
		#---? unclear argument type
		elif len(numpy.shape(args)) == 3 and numpy.shape(args)[0] == 1: args = tuple(args[0])
		elif numpy.shape(args)[1] > 2:
			#---if the arguments are just a list, it is assume that there are pairs of struct and slice
			args = tuple([[args[0][j],args[0][j+1]] for j in range(0,len(args[0]),2)])
		#---identify the keydims along which we will slice and what the slice should be
		for spec in args:
			keydims.append(self.struct[spec[0]])
			structs.append(spec[0])
			slices.append(spec[1])
		#---Generate a slice variable which will perform the final slice
		slicer = []
		#---construct the slice variable by iterating over each dimension and adding its particular slice
		for i in range(len(numpy.shape(self.data))):
			if i in keydims:
				if type(slices[keydims.index(i)]) != str: slicer.append(slices[keydims.index(i)])
				else: slicer.append(self.struct_opts[structs[keydims.index(i)]][slices[keydims.index(i)]])
			else: slicer.append(slice(0,None))
		#---the slicer object must be converted to a tuple before it is passed to the numpy array indexer
		slicer = tuple(slicer)
		#---return the data
		if flat: return flatten(numpy.array(self.data)[slicer])
		else: return numpy.array(self.data)[slicer]


#!/usr/bin/env python

import numpy
import collections
import itertools

#---master dictionary which provides specifications for different data types
master_datatypes = {
	'gr2d':{
		'struct':{'frame':0,'monolayer':1,'lipid':2},
		'struct_opts':None,
		'description':'Basic object which actually contains 3D radial distribution function (i.e. g(r)) data for AAMD bilayers. The root data is a normalized histogram or probability curve of the likelihood of reaching another particle. Particle specs are in the pkl name and the notes.',
		},
	'cells':{
		'struct':{'frame':0,'monolayer':1,'type':2},
		'struct_opts':{'type':{'points':0,'voronoi':1,'areas':2}},
		'description':'',
		},
	}

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
# NEEDS COMMENTS
	#---retrieve data
	def get(self,*args,**kwargs):
		print 'debug: get'
		#---Flatten the data as necessary, or depending on the calculation
		if 'flat' in kwargs.keys():
			print 'debug: flat'
			flat = kwargs['flat']
		else:
			print 'debug: flat else'
			print self.calctype
			if (self.calctype == 'grvoronoi' or self.calctype == 'gr'):
				flat = True
			elif (self.calctype == 'triangles' or self.calctype == 'cells'):
				flat = False
			else:
				flat = False
		print 'debug: flat ='
		print flat
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
		print 'slicer: '
		print slicer
		#---Return the data, flattened if desired
		if flat:
			return flatten(numpy.array(self.data)[slicer])
		else:
			print 'dbug: else return'
			print type(numpy.array(self.data)[slicer])
			print type(slicer)
			print type(self.data)
			tmpvar = numpy.array(self.data)[slicer]
			print tmpvar[0]
			tmpvar = numpy.array(self.data)
			print tmpvar[0]
			print self.data[0]
			print self.data[0][0]
			return numpy.array(self.data)[slicer]


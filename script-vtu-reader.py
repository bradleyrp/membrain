#!/usr/bin/python -i

#---this XML reader can import paraview-style VTU files 

from membrainrunner import *

location = ''
execfile('locations.py')

import numpy
import glob
import sys
import xml.etree.ElementTree as ET
import scipy
import os

need_trilist = 0

fname = '/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0'+\
	'/equilibrate/equib-conf-0.vtu'

#---Nb code from Ram via RWT to get the points
#---Nb RPB adapted the code to retrieve the induced_cur data
itemno = [j[1][1] for j in [i.items() for i in root[0].find('Piece').find('PointData').
	getchildren()]].index('induced_cur')
tree = ET.parse(fname)
root = tree.getroot()
for d in root.iter('UnstructuredGrid'):
	coord = d.find('Piece').find('Points').find('DataArray').text.split()
	coord = map(float,coord)
	n = coord.__len__()/3
	vcoord = (numpy.array(coord).reshape(n,3))
	c0s = d.find('Piece').find('PointData').getchildren()[itemno].text.split()
	c0s = map(float,c0s)
	if need_trilist:
		triangle = d.find('Piece').find('Cells').find('DataArray').text.split()
		triangle = map(int,triangle)
		n1 = triangle.__len__()/3
		trilist = (numpy.array(triangle).reshape(n1,3))
		data = [vcoord,trilist]


	


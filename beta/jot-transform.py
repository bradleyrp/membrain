#!/usr/bin/python -i

from membrainrunner import *
import numpy.linalg
execfile('plotter.py')

if 0:
	sel1 = array([[0,0,1],[1,0,1],[1,1,1]])
	#sel2 = array([[0,0,1],[0,1,2],[1,1,2]])
	meshpoints(transpose(sel1),color=(0,1,0),opacity=0.2)
	meshpoints(transpose(sel2),opacity=0.2,color=(1,0,0))
	translate = mean(sel1,axis=0)-mean(sel2,axis=0)
	meshpoints(transpose(sel1-translate),color=(0,0,1),opacity=0.2)
	#u, s, vh = numpy.linalg.svd(numpy.dot(sel2, sel1.T))
	#R = numpy.dot(u, vh)
	#sel3 = sel1 * R
	diff = sel2 - (sel1-translate)
	A = numpy.linalg.eig(diff)
	sel3 = sel1 * A[1]
	meshpoints(transpose(sel3), color=(0.7,0.7,0.7), opacity=1.0, scale_factor=1.5)
if 0:
	#---Following function translates and rotates three points using SVD.
	sel1 = array([[0,0,1],[1,0,1],[1,1,1]])
	sel2 = array([[0,0,2],[-1,0,2],[-1,-1,2]])
	#meshpoints(transpose(sel1),color=(0,1,0),opacity=0.2)
	#meshpoints(transpose(sel2),opacity=0.2,color=(1,0,0))
	#translate = mean(sel1,axis=0)-mean(sel2,axis=0)
	#sel3 = sel1-translate
	#meshpoints(transpose(sel3),color=(0,0,1),opacity=0.2)
	sel1b = sel1-mean(sel1,axis=0)
	sel2b = sel2-mean(sel2,axis=0)	
	A=dot(sel1b.T,sel2b)
	svdans = linalg.svd(A)
	m=numpy.diag((1,1,det(dot(svdans[2].T,svdans[0].T))))
	tmat=dot(dot((svdans[2].T),m),svdans[0].T)
	sel4 = array([list(dot(tmat,i)) for i in sel1b])
	meshpoints(transpose(sel1b),opacity=0.2,color=(1,0,0))
	meshpoints(transpose(sel2b),opacity=0.2,color=(0,1,0))
	meshpoints(transpose(sel4),color=(1,1,1),opacity=1)

'''
steps
read in a PI2P trajectory
identify the PI2P
take a phosphate group
move it
figure out how to save the new coordinates?
'''



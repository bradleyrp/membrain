#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')
#import numpy as N
#import pylab
#from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq
#import os

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

'''
TILEFILTER PROGRAM
Discretizes a membrane.
Options:
protein_subset_slice : slice of the stored protein points to use in the filter
height_direction : stores 1 or -1 as net-positive or net-negative tiles
cutoff_distance : when measuring around protein neighborhoods, how many grid-lengths to include
Notes:
Can specify different distance metrics in scipy.spatial.distance.cdist
'''


#---Analysis parameters
skip = 1
framecount = None
location = 'light'
execfile('locations-rpb.py')

#---Load
systemprefix_in = 'membrane-v614.md.part0002.skip10'
systemprefix = systemprefix_in+'.standard'
startpickle = pickles+'pkl.avgstruct.'+systemprefix_in+'.pkl'
heightfilter = -1
domainw = 4
manhatdist = 1

mset = unpickle(startpickle)
vecs = np.mean(mset.vecs,axis=0)

def g2d(params,x,y):
	'''Fitting function.'''
	z0,c0,x0,y0,sx,sy,th = params
	#return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	return z0+c0*exp(-(((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./w1**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	
def g2dresid(params,x,y,z,center):
	'''Fitting function residual.'''
	#return ((g2d(params,x,y)-z)**2)*(1 if center == None else 1./norm([x-center[0],y-center[1]]))
	#return ((g2d(params,x,y)-z)**2)/norm([x-center[0],y-center[1]])**2
	return ((g2d(params,x,y)-z)**2)
	
execfile('./jot-code/jot-sympy.py')

if 0:
	curvsmax = []
	params = []
	for fr in range(len(mset.surf)):
		print 'frame '+str(fr)
		#fr = 101
		cutoff_distance = 15.
		protein_subset_slice = slice(0,158)
		zvals={-1:0,1:1,0:2,2:3}
		which_brewer_colors = [1,5,0,4]
		colorcodes = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
		cutoff = cutoff_distance*10/(vecs[0]/mset.griddims[0])
		proteins = array(mset.protein)[:,protein_subset_slice]

		surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] 
			for i in range(mset.griddims[0]-1)])

		protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
		gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
		distmat = scipy.spatial.distance.cdist(protpos,gridpos)
		bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
			for i in range(len(distmat))] for w in v]))]
		buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
			dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
	
		target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete+buf==2)).T])
		pos = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete==1)).T])
		if 0:
			#meshplot(target,vecs=vecs)
			#meshplot(pos,vecs=vecs)
			#meshpoints(buf,scale_factor=10,vecs=vecs)
			meshplot(mset.surf[fr],vecs=vecs,show='surf',opacity=0.6)
			#meshpoints(target,vecs=vecs,scale_factor=20,color=(1,1,1))
			meshpoints(proteins[0]-[0,0,mean(mset.surf_position)],scale_factor=10,color=(1,0.24,0.59))
		target_com = [mean(target[:,0]),mean(target[:,1])]
		p_opt = leastsq(g2dresid,array([0,1,target_com[0],target_com[0],50,50,0]),args=(target[:,0],target[:,1],target[:,2],target_com))
		params.append(p_opt[0])
		fit = [[i[0],i[1],g2d(p_opt[0],i[0],i[1])] for i in target[:,0:2]]
		resid = [[i[0],i[1],g2d(p_opt[0],i[0],i[1])-i[2]] for i in target]
		resids1d = 10*(fit-target)[:,2]
		curvs = [10**4*gauss2dh(p_opt[0],i[0],i[1]) for i in target] ######### fix -1 problem
		if 0:
			meshpoints(fit,vecs=vecs,color=(1,1,1),scale_factor=resids1d)
			#meshpoints(target+[0,0,0],vecs=vecs,color=(1,1,1),scale_factor=curvs)
		curvs2 = [10*-1*gauss2dh(p_opt[0],i[0],i[1]) for i in target]
		#print max(curvs2)
		curvsmax.append(max(curvs2))
		#plt.hist([gauss2dh(p_opt[0],i[0],i[1]) for i in target]);plt.show()
	
if 1:
	validcurvs = [i for i in curvsmax if (i > 0.001 and i < 0.5)]
	validcurvsi = [i for i in range(len(curvsmax)) if (curvsmax[i] > 0.001 and curvsmax[i] < 0.5)]
	print 'log-averaged curvature is '+str(exp(mean(log(validcurvs))))
	fig = plt.figure()
	if 0:
		ax = fig.add_subplot(111)
		#ax.set_xscale('log', basey=10)
		ax.hist([log10(i) for i in validcurvs]);plt.show()
		plt.show()
	
	ax = fig.add_subplot(111)
	#ax.set_xscale('log', basey=10)
	ax.hist([log10(params[i,4]) for i in validcurvsi]);plt.show()
	plt.show()



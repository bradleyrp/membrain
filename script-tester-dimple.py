#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')

import os
from scipy.optimize import leastsq
import matplotlib as mpl
from pylab import *

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

griddims = [33,33]
vecs = array([ 655.04077148,  655.04077148,  249.67437744], dtype=float32)
startpts = [0,1,vecs[0]/2.,vecs[1]/2.,50,50,0]

execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	#return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2d_residual(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual.'''
	#return ((g2d(params,x,y)-z)**2)*(1 if center == None else 1./norm([x-center[0],y-center[1]]))
	#return ((g2d(params,x,y)-z)**2)/norm([x-center[0],y-center[1]])**2
	return ((gauss2d(params,x,y)-z)**2)

def test(params,show=False):
	dummysurf = array([[i*vecs[0]/(griddims[0]-1),j*vecs[1]/(griddims[1]-1),
		gauss2d(params,i*vecs[0]/(griddims[0]-1),
		j*vecs[1]/(griddims[1]-1))] 
		for j in range(griddims[1]-1) 
		for i in range(griddims[0]-1)])
	target=dummysurf
	target_com = [mean(target[:,0]),mean(target[:,1])]
	initpts = startpts
	initpts[2] = target_com[0]
	initpts[3] = target_com[1]
	p_opt = leastsq(gauss2d_residual,array(initpts),args=(target[:,0],target[:,1],target[:,2]))
	print list([round(i,3) for i in p_opt[0]])
	print list([round(i,3) for i in params])
	if show:
		fig = plt.figure(figsize=(8,10))
		zmap = array([[gauss2d(params,i*vecs[0]/(griddims[0]-1),
			j*vecs[1]/(griddims[1]-1)) 
			for j in range(griddims[1]-1)]
			for i in range(griddims[0]-1)])
		ax = plt.subplot(211)
		img = ax.imshow(zmap,interpolation='nearest',origin='LowerLeft',)
		cax = inset_axes(ax,
			width="5%",
			height="100%",
			bbox_transform=ax.transAxes,
			bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
			loc= 1)
		fig.colorbar(img,cax=cax)
		hmap = array([[-1*gauss2dh(params,i*vecs[0]/(griddims[0]-1),
			j*vecs[1]/(griddims[1]-1))
			for j in range(griddims[1]-1)]
			for i in range(griddims[0]-1)])
		ax2 = plt.subplot(212)
		img2 = ax2.imshow(hmap,interpolation='nearest',origin='LowerLeft',)
		cax2 = inset_axes(ax2,
			width="5%",
			height="100%",
			bbox_transform=ax2.transAxes,
			bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
			loc=1)
		fig.colorbar(img2,cax=cax2)
		plt.show()

#---trials
paramtrials = [
	[0,100.,vecs[0]/2.,vecs[1]/2.,100,200,0],
	[0,-100.,vecs[0]/2.,vecs[1]/2.,100,200,0],
	[0,-100.,vecs[0]/2.,vecs[1]/2.,100,200,pi/4.],
	[0,10.,vecs[0]/2.,vecs[1]/2.,10,20,pi/4.],
	[0,10.,vecs[0]/2.,vecs[1]/2.,20,40,pi/4.],
	[0,10.,vecs[0]/2.,vecs[1]/2.,40,80,pi/4.],
	[0,-10.,vecs[0]/2.,vecs[1]/2.,40,80,pi/4.],
	[0,100.,vecs[0]/2.,vecs[1]/2.,100,200,0],
	[0,200.,vecs[0]/2.,vecs[1]/2.,100,200,pi/4.]]
test(paramtrials[-1],show=True)


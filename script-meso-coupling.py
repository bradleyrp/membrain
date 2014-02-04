#!/usr/bin/python -i

if 0:

	from membrainrunner import *

	import numpy
	import glob
	import sys
	import xml.etree.ElementTree as ET
	import scipy
	import os

	#---SETTINGS
	#-------------------------------------------------------------------------------------------------------------

	#---methods
	location = ''
	execfile('locations.py')

	#---settings
	vtudir = '/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/'
	xyzform = 'rect'
	nbase = 22
	length = None

	#---FUNCTIONS
	#-------------------------------------------------------------------------------------------------------------

	def fftredundant(dat):
		# was divided by lenscale before?
		return fft.fftshift(fft.fft2(array(dat)[:-1,:-1]))

	#---MAIN
	#-------------------------------------------------------------------------------------------------------------

	extradat = mset.load_points_vtu(vtudir,extra_props='induced_cur',start=100,end=110,nbase=nbase)
	mset.surfacer()
	c0s = mset.surfacer_general(array(extradat)[:,0])
	# plot with: meshplot(array([[mset.xyzs[0][i][0],mset.xyzs[0][i][1],extradat[0][i]] for i in range(574)]))

if 0:
	grid = mset.griddims
	qs0 = [[x,y] for x in linspace(0,mset.vecs[0][0],mset.griddims[0]) 
		for y in linspace(0,mset.vecs[0][1],mset.griddims[1])]
	qs1 = [[x,y] for x in linspace(0,mset.vecs[0][0],mset.griddims[0]) 
		for y in linspace(0,mset.vecs[0][1],mset.griddims[1])]
	qis0 = [[x,y] for x in range(1,grid[0]-1) 
		for y in range(1,grid[1]-1)]
	qis1 = [[x,y] for x in range(1,grid[0]-1) 
		for y in range(1,grid[1]-1)]
	hq = fftredundant(mset.surf[0])
	cq = fftredundant(c0s[0])
	dat_raw = [[arccos(np.dot(q0,q1)/np.linalg.norm(q0)/np.linalg.norm(q1)),cq[q0[0],q0[1]]*hq[q1[0],q1[1]]] for q0 in qis0[1:] for q1 in qis1[1:]]
	dat = [[float(abs(i[0])),float(abs(i[1]))] for i in dat_raw if (not np.isnan(abs(i[0])) and not np.isnan(abs(i[1])))]
	rangey=(min(array(dat)[:,1]),max(array(dat)[:,1]))
	rangex=(min(array(dat)[:,0]),max(array(dat)[:,0]))
if 1:
	import matplotlib.gridspec as gridspec
	fig = plt.figure(figsize=(6,6))
	gs = gridspec.GridSpec(1,1)
	ax = plt.subplot(gs[0])
	fingerprint, xedges, yedges = numpy.histogram2d(array(dat)[:,0],array(dat)[:,1],bins=20,range=([rangex[0],rangex[1]],[rangey[0],10**-1]))
	midx = (xedges[1:]+xedges[:-1])/2.
	midy = (yedges[1:]+yedges[:-1])/2.
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	#cmap.set_bad(cmap(0),1.)
	ax.imshow(array(fingerprint).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=None,cmap=cmap)
	plt.show()
	#[[q0,q1] for q0 in qis0 for q1 in qis1]
	
			




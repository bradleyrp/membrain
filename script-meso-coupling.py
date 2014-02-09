#!/usr/bin/python -i

from membrainrunner import *

import numpy
import glob
import sys
import xml.etree.ElementTree as ET
import scipy
import os
from numpy.linalg import norm
import matplotlib.gridspec as gridspec

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---methods
location = ''
execfile('locations.py')

#---settings
xyzform = 'rect'
nbase = 22
length = None

#---plan
analysis_descriptors = [
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-anis-22-run1-rep-0-1500-2000',True,''),
	('/home/rpb/worker/repo-membrane/membrane-v2002/bare-rep-0/equilibrate',1500,2000,
		'v2002-bare-rep-0-1500-2000',True,''),
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl')]
ad = analysis_descriptors[-1]
vtudir,start,end,testname,ismeso,pklname = ad

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def fftredundant(dat):
	return fft.fftshift(fft.fft2(array(dat)[:-1,:-1]))

#---MAIN
#-------------------------------------------------------------------------------------------------------------

	#|------load
	#||-----calc
	#|||----calc2
	#||||---plot
seq='0001'
	   	
if int(seq[0]):
	if ismeso:
		extradat = mset.load_points_vtu(vtudir,extra_props='induced_cur',start=start,end=end,nbase=nbase)
		mset.surfacer()
		c0s = mset.surfacer_general(array(extradat)[:,0])
	else:
		mset = unpickle(pickles+pklname)
if int(seq[1]):
	grid = mset.griddims
	qxn,qyn = grid
	q0x,q0y = [int(round(i/2.))-1 for i in grid]
	q = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qp = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	q0x,q0y = [int(round(i/2.))-1 for i in shape(q)[0:2]]
	
	hqs = []
	for i in range(len(mset.surf)):
		hqs.append(fftredundant(mset.surf[i]))
	cq = fftredundant(c0s[0])
	print 0
	term0 = [[
		norm(sum([norm(q[qjx][qjy])**2*norm(qp[qix][qiy])**2*mean([hqs[i][qjx][qjy]*hqs[i][qix][qiy] 
		for i in range(len(hqs))])
		for qix in range(qxn-1) 
		for qiy in range(qyn-1)]))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	print 1
	term1 = [[
		norm(sum([norm(q[qjx][qjy])**2*mean([hqs[i][qjx][qjy]*cq[qix][qiy] for i in range(len(hqs))])
		for qix in range(qxn-1) 
		for qiy in range(qyn-1)]))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	print 2
	term2 = [[
		norm(sum([norm(qp[qix][qiy])**2*mean([cq[qjx][qjy]*hqs[i][qix][qiy] for i in range(len(hqs))])
		for qix in range(qxn-1) 
		for qiy in range(qyn-1)]))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	print 3
	term3 = [[
		norm(sum([mean([cq[qjx][qjy]*cq[qix][qiy] for i in range(len(hqs))])
		for qix in range(qxn-1) 
		for qiy in range(qyn-1)]))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
if int(seq[2]):
	grid = mset.griddims
	qxn,qyn = grid
	q0x,q0y = [int(round(i/2.))-1 for i in grid]
	q = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qp = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	q0x,q0y = [int(round(i/2.))-1 for i in shape(q)[0:2]]
	
	hqs = []
	for i in range(len(mset.surf)):
		hqs.append(fftredundant(mset.surf[i]))
	cq = fftredundant(c0s[0])
	print 0
	term0 = [[
		norm(sum([norm(q[qjx][qjy])**2*norm(qp[qix][qiy])**2*mean(np.real([hqs[i][qjx][qjy]*hqs[i][qix][qiy]
		for i in range(len(hqs))]))
		for qix in range(qxn-1) 
		for qiy in range(qyn-1)]))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	if ismeso:
		print 1
		term1 = [[
			norm(sum([norm(q[qjx][qjy])**2*mean(np.real([hqs[i][qjx][qjy]*cq[qix][qiy] for i in range(len(hqs))]))
			for qix in range(qxn-1)
			for qiy in range(qyn-1)]))
			for qjy in range(qyn-1)]
			for qjx in range(qxn-1)]
		print 2
		term2 = [[
			norm(sum([norm(qp[qix][qiy])**2*mean(np.real([cq[qjx][qjy]*hqs[i][qix][qiy] for i in range(len(hqs))]))
			for qix in range(qxn-1)
			for qiy in range(qyn-1)]))
			for qjy in range(qyn-1)]
			for qjx in range(qxn-1)]
		print 3
		term3 = [[
			norm(sum([mean(np.real([cq[qjx][qjy]*cq[qix][qiy] for i in range(len(hqs))]))
			for qix in range(qxn-1) 
			for qiy in range(qyn-1)]))
			for qjy in range(qyn-1)]
			for qjx in range(qxn-1)]
	else:
		term1 = term0
		term2 = term0
		term3 = term0
if int(seq[3]):
	islognorm = True
	fig = plt.figure(figsize=(20,12))
	cmap = mpl.cm.jet
	gs = gridspec.GridSpec(2,4)
	ax = plt.subplot(gs[0,0])
	ax.set_title('term 1')
	ax.imshow(array(term0).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term0).max()))
	ax = plt.subplot(gs[0,1])
	ax.set_title('term 2')
	ax.imshow(array(term1).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term1).max()))
	ax = plt.subplot(gs[1,0])
	ax.set_title('term 3')
	ax.imshow(array(term2).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term2).max()))
	ax = plt.subplot(gs[1,1])
	ax.set_title('term 4')
	ax.imshow(array(term3).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term3).max()))
	
	ax = plt.subplot(gs[0,2])
	ax.set_title('sum')
	ax.imshow(array(array(term0)-array(term1)-array(term2)+array(term3)).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(array(term0)-array(term1)-array(term2)+array(term3)).max()))

	ax = plt.subplot(gs[1,2])
	ax.set_title('sum raw ')
	ax.imshow(array(array(term0)-array(term1)-array(term2)+array(term3)).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap)
	ax.set_xlabel(str(array(array(term0)-array(term1)-array(term2)+array(term3)).max()))
	plt.savefig(pickles+'fig-bilayer-couple-view'+testname+('-lognorm' if islognorm else '')+'.png',
		dpi=500,bbox_inches='tight')	
	plt.show()

			




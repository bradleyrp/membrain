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
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-bare-22-run1-rep-0-1500-2000',True,''),
	('/home/rpb/worker/repo-membrane/membrane-v2002/bare-rep-0/equilibrate',1500,2000,
		'v2002-bare-rep-0-1500-2000',True,''),
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl')]
ad = analysis_descriptors[1]
vtudir,start,end,testname,ismeso,pklname = ad

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def fftredundant(dat):
	return fft.fftshift(fft.fft2(array(dat)[:-1,:-1]))
	#return fft.fft2(array(dat)[:-1,:-1])

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#------|------load
#------||-----calc-v1
#------|||----calc-v2
#------||||---plot
seq = '0011'
	   	
if int(seq[0]):
	if ismeso:
		extradat = mset.load_points_vtu(vtudir,extra_props='induced_cur',start=start,end=end,nbase=nbase)
		mset.surfacer()
		c0s = mset.surfacer_general(array(extradat)[:,0])
	else:
		mset = unpickle(pickles+pklname)
#---METHOD A. first method, probably incorrectly includes sum over q'
if int(seq[1]):
	grid = mset.griddims
	qxn,qyn = grid
	q0x,q0y = [int(round(i/2.))-1 for i in grid]
	q0x,q0y = [float(q0x)/grid[0],float(q0y)/grid[1]]*array(mean(mset.vecs,axis=0))[0:2]
	q = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qp = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	hqs = []
	for i in range(len(mset.surf)):
		hqs.append(fftredundant(mset.surf[i]))
	if not ismeso:
		c0s = [zeros(grid)]
	cq = fftredundant(zeros(grid))
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
			norm(sum([norm(q[qjx][qjy])**2*mean(np.real([hqs[i][qjx][qjy]*cq[qix][qiy] 
			for i in range(len(hqs))]))
			for qix in range(qxn-1)
			for qiy in range(qyn-1)]))
			for qjy in range(qyn-1)]
			for qjx in range(qxn-1)]
		print 2
		term2 = [[
			norm(sum([norm(qp[qix][qiy])**2*mean(np.real([cq[qjx][qjy]*hqs[i][qix][qiy] 
			for i in range(len(hqs))]))
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
#---METHOD B. second method collapsed q
if int(seq[2]):
	mset.calculate_undulation_spectrum()
	mset.analyze_undulations()
	grid = mset.griddims
	qxn,qyn = grid
	q0x,q0y = [int(round(i/2.))-0 for i in grid]
	q0x,q0y = [float(q0x)/grid[0],float(q0y)/grid[1]]*array(mean(mset.vecs,axis=0))[0:2]
	#q0x,q0y = 0,0
	q = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qp = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qmags = [[norm([x-mset.vecs[0][0]/2,y-mset.vecs[0][1]/2]) for y in linspace(0,mset.vecs[0][1],grid[1])] 
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qmags[12][14] = 0.
	hqs = []
	for i in range(len(mset.surf)):
		hqs.append((fftredundant(mset.surf[i])/double(grid[0]*grid[1]))**2)
	if not ismeso:
		c0s = [zeros(grid)]
	cq = (fftredundant(c0s[0])/double(grid[0]*grid[1]))**2
	print 0
	term0 = [[
		qmags[qjx][qjy]**4*norm(mean(np.real([hqs[i][qjx][qjy]*hqs[i][qjx][qjy]
		for i in range(len(hqs))])))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	print 1
	term1 = [[
		qmags[qjx][qjy]**2*norm(mean(np.real([cq[qjx][qjy]*hqs[i][qjx][qjy]
		for i in range(len(hqs))])))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]	
	print 2
	term2 = [[
		qmags[qjx][qjy]**2*norm(mean(np.real([hqs[i][qjx][qjy]*cq[qjx][qjy]
		for i in range(len(hqs))])))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
	print 3
	term3 = [[
		norm(mean(np.real([cq[qjx][qjy]*cq[qjx][qjy]
		for i in range(len(hqs))])))
		for qjy in range(qyn-1)]
		for qjx in range(qxn-1)]
if int(seq[3]):
	islognorm = True
	fig = plt.figure(figsize=(18,8))
	cmap = mpl.cm.jet
	cmap = mpl.cm.get_cmap('RdGy',100)
	cmap = mpl.cm.binary
	gs = gridspec.GridSpec(2,5)
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
	
	termsum = array(term0)-array(term1)-array(term2)+array(term3)

	ax = plt.subplot(gs[0,2])
	ax.set_title('sum')
	ax.imshow(array(termsum).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(array(term0)-array(term1)-array(term2)+array(term3)).max()))

	ax = plt.subplot(gs[1,2])
	ax.set_title('sum raw ')
	ax.imshow(array(termsum).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap)
	ax.set_xlabel(str(array(array(term0)-array(term1)-array(term2)+array(term3)).max()))
	
	
	ax2 = plt.subplot(gs[:,3:])
	
	lenscale = 1.0
	qmagfilter=[10**-6,10**6]
	#ax2.set_xlabel(r"$\left|q\right|$",fontsize=18)
	#ax2.set_ylabel(r"$\left\langle z_{q}z_{-q}\right\rangle$",fontsize=18)
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.grid(True)
	spectrum = array(sorted(zip(*[mset.qrawmean,mset.uqrawmean,mset.uqrawstd]), key=lambda x: x[0]))[1:]
	ax2.scatter(spectrum[:,0],spectrum[:,1],marker='o',color='k')
	spectrumf = array(filter(lambda x: x[0] >= qmagfilter[0] and x[0] <= qmagfilter[1], spectrum))
	spectrumf2 = array(filter(lambda x: x[0] >= 1./10.*qmagfilter[0] and x[0] <= qmagfilter[1]*10., spectrum))
	[bz,az]=numpy.polyfit(log(spectrumf[:,0]),log(spectrumf[:,1]),1)
	area = double(mean([i[0]*i[1] for i in mset.vecs])/lenscale**2)
	print 'q-magnitude scaling = '+str(bz)
	kappa = 1/exp(az)/area
	print 'kappa = '+str(kappa)
	leftcom = [mean(log(spectrumf[:,0])),mean(log(spectrumf[:,1]))]
	az_enforced = leftcom[1]+4.*leftcom[0]
	kappa_enforced = 1./exp(az_enforced)/area
	print 'kappa_enforced = '+str(kappa_enforced)
	ymod=[exp(az)*(i**bz) for i in spectrumf[:,0]]
	xmod=[i for i in spectrumf[:,0]]
	ymod2=[exp(az)*(i**bz) for i in spectrumf2[:,0]]
	xmod2=[i for i in spectrumf2[:,0]]
	ymod3=[exp(az)*(i**bz) for i in spectrumf2[:,0]]
	xmod3=[i for i in spectrumf2[:,0]]
	ymod4=[exp(az_enforced)*(i**-4.) for i in spectrumf2[:,0]]
	xmod4=[i for i in spectrumf2[:,0]]
	ax2.plot(xmod2,ymod2,color='#FF3399',linestyle='-',linewidth=2.5)
	ax2.plot(xmod4,ymod4,color='#3399FF',linestyle='-',linewidth=2.5)
	ax2.text(0.1, 0.2, r'$\kappa = %3.2f$'%kappa, transform=ax2.transAxes,fontsize=16)
	
	unscaled0 = array(term0)/array(qmags)[:-1,:-1]
	specs0 = array([[qmags[qjx][qjy],unscaled0[qjx][qjy]/qmags[qjx][qjy]**4]
		for qjy in range(qyn-1)
		for qjx in range(qxn-1) if (qjx != 12 and qjy != 14)])
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.scatter(specs0[:,0],specs0[:,1],color='r')
	unscaled1 = array(term1)/array(qmags)[:-1,:-1]
	specs1 = array([[qmags[qjx][qjy],unscaled1[qjx][qjy]/qmags[qjx][qjy]**2]
		for qjy in range(qyn-1)
		for qjx in range(qxn-1) if (qjx != 12 and qjy != 14)])
	ax2.scatter(specs1[:,0],specs1[:,1],color='b')
	unscaled3 = array(term3)/array(qmags)[:-1,:-1]
	specs3 = array([[qmags[qjx][qjy],unscaled3[qjx][qjy]]
		for qjy in range(qyn-1)
		for qjx in range(qxn-1) if (qjx != 12 and qjy != 14)])
	ax2.scatter(specs3[:,0],specs3[:,1],color='g')

	ax2.set_ylim((0.9*array(term3).min(),1.1*max(array(term0).max(),max(spectrum[:,1]))))
	ax2.set_xlim((0.9*min(spectrum[:,0]),1.1*sort(array(qmags).max())))
	


	plt.savefig(pickles+'fig-bilayer-couple-view-'+testname+('-lognorm' if islognorm else '')+'.png',
		dpi=500,bbox_inches='tight')
	
	plt.show()


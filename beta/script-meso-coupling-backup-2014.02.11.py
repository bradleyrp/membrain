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
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl')]
ad = analysis_descriptors[0]
vtudir,start,end,testname,ismeso,pklname = ad

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def fftredundant(dat):
	return fft.fftshift(fft.fft2(array(dat)[:-1,:-1]))
	#return fft.fft2(array(dat)[:-1,:-1])

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#------|------load
#------||-----calc
#------|||----test, 2d
#------||||---plot
seq = '0111'

if int(seq[0]) or mset.surf == []:
	if ismeso:
		extradat = mset.load_points_vtu(vtudir,extra_props='induced_cur',start=start,end=end,nbase=nbase)
		mset.surfacer()
		c0s = mset.surfacer_general(array(extradat)[:,0])
	else:
		mset = unpickle(pickles+pklname)
#---calculations, trial B
if int(seq[1]):
	mset.calculate_undulation_spectrum(redundant=1)
	mset.analyze_undulations()
	grid = mset.griddims
	qxn,qyn = grid
	q0x,q0y = [int(round(i/2.))-1 for i in grid]
	#q0x,q0y = [float(q0x)/grid[0],float(q0y)/grid[1]]*array(mean(mset.vecs,axis=0))[0:2]
	Lx,Ly = mean(mset.vecs,axis=0)[0:2]
	#q0x,q0y = 0,0
	q = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qp = [[[x-q0x,y-q0y] for y in linspace(0,mset.vecs[0][1],grid[1])]
		for x in linspace(0,mset.vecs[0][0],grid[0])]
	qmags = [[norm([x-mset.vecs[0][0]/2,y-mset.vecs[0][1]/2]) for y in linspace(0,mset.vecs[0][1],grid[1])] for x in linspace(0,mset.vecs[0][0],grid[0])]
	#q0x,q0y = Lx*(grid[0]-0)/(grid[0]+1)/2.,Ly*(grid[1]-0)/(grid[1]+1)/2.
	qmags = array([[sqrt(((i-q0x)/((Lx)/1.)*2*pi)**2+((j-q0y)/((Ly)/1.)*2*pi)**2) for j in range(0,grid[1])] for i in range(0,grid[0])])
	#qmags[12][14] = 0.
	hqs = []
	for i in range(len(mset.surf)):
		hqs.append(fftredundant(mset.surf[i]))
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
#---check 2D FFTs
if int(seq[2]):
	redundant = 1
	m,n = grid[0]-1,grid[1]-1
	center = [12,14]

	uqrawmean = mean((abs(array(hqs))/double(m*n))**2,axis=0)
	qraw = qmags
	spectrum1d = array([[qmags[i,j],uqrawmean[i,j]] for j in range(n) for i in range(m)])
	
	#---ATOMIC OPERATION !!!
	hqsa = array(hqs)[:,-1:center[0]%2:-1,-1:center[0]%2:-1]
	hqsb = array(hqs)[:,(center[0]+1)%2:,(center[0]%2+1):]
	uqrawmeanlit = mean(abs(hqsa*hqsb)/double(m*n)**2,axis=0)
	c = center[0]-(center[0]+1)%2,center[1]-(center[1]+1)%2
	spectrum1d = array([[qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)],uqrawmeanlit[i,j]] for j in range(n-((center[1]+1)%2)) for i in range(m-((center[0]+1)%2))])
	
	cqb = cq[(center[0]+1)%2:,(center[0]%2+1):]
	cqa = cq[-1:center[0]%2:-1,-1:center[0]%2:-1]
	term1alt = array([[norm(mean(np.real([cqb[qjx][qjy]*hqsa[i][qjx][qjy] for i in range(len(hqs))]))) for qjy in range(qyn-1-(center[0]+1)%2)] for qjx in range(qxn-1-(center[0]+1)%2)])
	term1spec = array([[qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)],term1alt[i,j]] for j in range(n-((center[1]+1)%2)) for i in range(m-((center[0]+1)%2))]) 
	term2alt = array([[norm(mean(np.real([hqsb[i][qjx][qjy]*cqa[qjx][qjy] for i in range(len(hqs))]))) for qjy in range(qyn-1-(center[0]+1)%2)] for qjx in range(qxn-1-(center[0]+1)%2)])
	term2spec = array([[qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)],term2alt[i,j]] for j in range(n-((center[1]+1)%2)) for i in range(m-((center[0]+1)%2))]) 
	term3alt = array([[norm(mean(np.real([cqb[qjx][qjy]*cqa[qjx][qjy] for i in range(len(hqs))]))) for qjy in range(qyn-1-(center[0]+1)%2)] for qjx in range(qxn-1-(center[0]+1)%2)])
	term3spec = array([[qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)],term3alt[i,j]] for j in range(n-((center[1]+1)%2)) for i in range(m-((center[0]+1)%2))]) 

	termsumalt2d = [[uqrawmeanlit[i,j]*qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)]**4-
		term1alt[i,j]*qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)]**2-
		term2alt[i,j]*qmags[i+((center[0]+1)%2),j+((center[1]+1)%2)]**2+
		term3alt[i,j] for j in range(n-1)] for i in range(m-1)]
	termsumalt = [spectrum1d[i,1]*spectrum1d[i,0]**4-2*term1spec[i,1]*term1spec[i,0]**2+term3spec[i,1] for i in range(len(spectrum1d))]

	
#---plot comparison
if int(seq[3]):
	islognorm = True
	fig = plt.figure(figsize=(18,8))
	cmap = mpl.cm.jet
	cmap = mpl.cm.get_cmap('RdGy',100)
	cmap = mpl.cm.binary
	gs = gridspec.GridSpec(2,5)
	ax = plt.subplot(gs[0,0])
	ax.set_title('term 1')
	ax.imshow(array(uqrawmeanlit).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(uqrawmeanlit).max()))
	ax = plt.subplot(gs[0,1])
	ax.set_title('term 2')
	ax.imshow(array(term1alt).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term1alt).max()))
	ax = plt.subplot(gs[1,0])
	ax.set_title('term 3')
	ax.imshow(array(term2alt).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term2alt).max()))
	ax = plt.subplot(gs[1,1])
	ax.set_title('term 4')
	ax.imshow(array(term3alt).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(term3alt).max()))
	
	termsum = array(term0)-array(term1)-array(term2)+array(term3)

	ax = plt.subplot(gs[0,2])
	ax.set_title('sum')
	ax.imshow(array(termsumalt2d).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(termsumalt2d).max()))

	ax = plt.subplot(gs[1,2])
	ax.set_title('sum raw ')
	ax.imshow(array(termsumalt2d).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap)
	ax.set_xlabel(str(array(termsumalt2d).max()))

	ax2 = plt.subplot(gs[:,3:])
	
	lenscale = 1.0
	qmagfilter=[10**-10,10**6]
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
	if 0:
		unscaled0 = sqrt(array(array(term0)))
		specs0 = array([[qmags[qjx][qjy],unscaled0[qjx][qjy]/qmags[qjx][qjy]**0] 
			for qjy in range(qyn-1) for qjx in range(qxn-1)])
		ax2.set_xscale('log')
		ax2.set_yscale('log')
		ax2.scatter(specs0[:,0],specs0[:,1],color='r')
		unscaled1 = sqrt(array(array(term1)))
		specs1 = array([[qmags[qjx][qjy],unscaled1[qjx][qjy]/qmags[qjx][qjy]**0]
			for qjy in range(qyn-1)
			for qjx in range(qxn-1)])
		ax2.scatter(specs1[:,0],specs1[:,1],color='b')
		unscaled3 = sqrt(array(array(term3)))
		specs3 = array([[qmags[qjx][qjy],unscaled3[qjx][qjy]]
			for qjy in range(qyn-1)
			for qjx in range(qxn-1)])
		ax2.scatter(specs3[:,0],specs3[:,1],color='g')
	
		unscaledall = sqrt(array(array(termsum)))
		specsall = array([[qmags[qjx][qjy],unscaledall[qjx][qjy]]
			for qjy in range(qyn-1)
			for qjx in range(qxn-1)])
		ax2.scatter(specsall[:,0],specsall[:,1],color='c')
		#ax2.set_ylim((0.9*array(term3).min(),1.1*max(array(term0).max(),max(spectrum[:,1]))))
		#ax2.set_xlim((0.9*min(spectrum[:,0]),1.1*sort(array(qmags).max())))
	else:
		ax2.scatter(spectrum1d[:,0],spectrum1d[:,1],color='c')
	#ax2.set_ylim((10**-14,10**1))
	#ax2.set_xlim((10**-1,10**1))
	ax2.scatter(term1spec[:,0],term1spec[:,1],color='r')	
	ax2.scatter(term2spec[:,0],term2spec[:,1],color='b')
	ax2.scatter(term3spec[:,0],term3spec[:,1],color='g')
	ax2.scatter(term3spec[:,0],termsumalt,color='m')
	plt.savefig(pickles+'fig-bilayer-couple-view-'+testname+('-lognorm' if islognorm else '')+'.png',
		dpi=500,bbox_inches='tight')
	plt.show()

#---EXTRAS
#-------------------------------------------------------------------------------------------------------------
	
#---calculation, trial A
if 0:
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



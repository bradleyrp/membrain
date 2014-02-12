#!/usr/bin/python -i

from membrainrunner import *

import numpy, glob, sys, scipy, os
import xml.etree.ElementTree as ET
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
		1500,1750,'v2002-t2-anis-22-run1-rep-0-1500-2000',True,''),
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,1510,'v2002-t2-bare-22-run1-rep-0-1500-2000',True,''),
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl')]
ad = analysis_descriptors[0]
vtudir,start,end,testname,ismeso,pklname = ad

#---method
#------|------load
#------||-----calc
#------|||----plot
seq = '011'

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def fftwrap(dat,redundant=1):
	'''This function wraps the standard discrete FFT for a system with possible-redundant rows.'''
	trim = -1 if redundant == 1 else None 
	return fft.fftshift(fft.fft2(array(dat)[:trim,:trim]))
	
def autocorr(dat,direct=0):
	'''Standard procedure for turning a possibly even grid into an odd one with a distinct center.'''
	m,n = shape(dat[0])
	if direct == 0:
		return array(dat)[:,slice((1 if m%2==0 else None),None),
			slice((1 if n%2==0 else None),None)]
	elif direct == 1:
		return array(dat)[:,slice((-1 if m%2==1 else None),(0 if m%2==0 else None),-1),
			slice((-1 if n%2==1 else None),(0 if n%2==0 else None),-1)]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load and interpolate
if int(seq[0]) or mset.surf == []:
	if ismeso:
		c0sraw = array(mset.load_points_vtu(vtudir,extra_props='induced_cur',
			start=start,end=end,nbase=nbase))[:,0]
		mset.surfacer()
		c0s = mset.surfacer_general(c0sraw)
	else:
		mset = unpickle(pickles+pklname)

#---calculate mode couplings
if int(seq[1]):
	mset.calculate_undulation_spectrum()
	mset.analyze_undulations()
	grid = mset.griddims
	m,n = grid[0]-1,grid[1]-1
	hqs = [fftwrap(mset.surf[i])/double(m*n) for i in range(len(mset.surf))]
	cqs = [fftwrap(c0s[i])/double(m*n) for i in range(len(c0s))]
	Lx,Ly = mean(mset.vecs,axis=0)[0:2]
	cm,cn = [int(round(i/2.-1)) for i in shape(hqs[0])]
	qmags = array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
		for j in range(0,n)] for i in range(0,m)])
	hqsa = autocorr(hqs,direct=0)
	hqsb = autocorr(hqs,direct=1)
	center = [cm,cn]
	cqsa = autocorr(cqs,direct=0)
	cqsb = autocorr(cqs,direct=1)
	qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]
	mt,nt = shape(hqsa[0])
	t0 = uqrawmean = mean((abs(array(hqsa)))**2,axis=0)
	spectrum1d = array([[qmagst[i,j],uqrawmean[i,j]] for j in range(nt) for i in range(mt)])
	t1 = mean(abs(hqsa)*abs(cqsb),axis=0)
	t2 = mean(abs(cqsa*hqsb),axis=0)
	t3 = mean(abs(cqsa*cqsb),axis=0)
	t1spec = array([[qmagst[i,j],t1[i,j]] for j in range(nt) for i in range(mt)])
	t2spec = array([[qmagst[i,j],t2[i,j]] for j in range(nt) for i in range(mt)])
	t3spec = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])
	termsum = t0*qmagst**4-t1*qmagst**2-t2*qmagst**2+t3
	termsumspec = array([[qmagst[i,j],termsum[i,j]] for j in range(nt) for i in range(mt)])
	t0spec2 = array([[qmagst[i,j],uqrawmean[i,j]*qmagst[i,j]**4] for j in range(nt) for i in range(mt)])

#---plots
if int(seq[2]):
	islognorm = True
	fig = plt.figure(figsize=(18,8))
	cmap = mpl.cm.jet
	gs = gridspec.GridSpec(2,5)
	ax = plt.subplot(gs[0,0])
	ax.set_title('term 1')
	ax.imshow(array(t0).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(t0).max()))
	ax = plt.subplot(gs[0,1])
	ax.set_title('term 2')
	ax.imshow(array(t1).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(t1).max()))
	ax = plt.subplot(gs[1,0])
	ax.set_title('term 3')
	ax.imshow(array(t2).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(t2).max()))
	ax = plt.subplot(gs[1,1])
	ax.set_title('term 4')
	ax.imshow(array(t3).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(t3).max()))
	
	#termsum = array(term0)-array(term1)-array(term2)+array(term3)

	ax = plt.subplot(gs[0,2])
	ax.set_title('sum')
	ax.imshow(array(termsum).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None))
	ax.set_xlabel(str(array(termsum).max()))

	ax = plt.subplot(gs[1,2])
	ax.set_title('sum raw ')
	ax.imshow(array(termsum).T, extent=None,
		interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap)
	ax.set_xlabel(str(array(termsum).max()))

	#---combined 1D spectra
	ax2 = plt.subplot(gs[:,3:])

	#---plot undulation spectra 1D
	lenscale = 1.0
	qmagfilter=[10**-10,10**6]
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
	ax2.scatter(spectrum1d[:,0],spectrum1d[:,1],color='c')
	ax2.scatter(t1spec[:,0],t1spec[:,1],color='b',marker='+')	
	ax2.scatter(t2spec[:,0],t2spec[:,1],color='p',marker='x')
	ax2.scatter(t3spec[:,0],t3spec[:,1],color='r')
	ax2.scatter(t0spec2[:,0],t0spec2[:,1],color='g')
	ax2.scatter(termsumspec[:,0],termsumspec[:,1],color='m')
	plt.savefig(pickles+'fig-bilayer-couple-view-'+testname+('-lognorm' if islognorm else '')+'.png',
		dpi=500,bbox_inches='tight')
	plt.show()
	

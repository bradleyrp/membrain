#!/usr/bin/python -i

from membrainrunner import *

import numpy, glob, sys, scipy, os
import xml.etree.ElementTree as ET
from numpy.linalg import norm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---methods
location = ''
execfile('locations.py')

#---settings
xyzform = 'rect'
nbase = 22
lenscale = 4.0
length = None
qmagfilter=[10**-10,10**6]

#---plan
analysis_descriptors = [
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t2-anis-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-anis-22-run1-rep-0-1500-2000',True,'',False),
	('/home/rpb/worker/repo-membrane/mesoscale-v2002/t1-bare-22/run1-size-sweep/rep-0/equilibrate/',
		1500,2000,'v2002-t2-bare-22-run1-rep-0-1500-2000',True,'',True),
	('',0,0,'v700',False,'pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',False)]
ad = analysis_descriptors[0]
vtudir,start,end,testname,ismeso,pklname,isbare = ad
barecompare = True
baresys = 1

#---method
    '''|---------load
       ||--------calc
       |||-------plot terms
       ||||------plot compare hqhq
       |||||-----plot compare hqhqq4
       ||||||----plot spectrum, 1D '''
seq = '010001'

#---settings
showplots = True
sskip=4
clrs = [(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)[i] for i in [0,2,3,4,1,5,6,7]]

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
			start=start,end=end,nbase=nbase,lenscale=lenscale))[:,0]
		mset.surfacer()
		c0s = mset.surfacer_general(c0sraw)
		#---note that you have to fix the scaling here. lenscale is used in rezipgrid, problem for c0
	else:
		mset = unpickle(pickles+pklname)
	if barecompare:
		vtudirb,startb,endb,testnameb,ismesob,pklnameb,isbareb = analysis_descriptors[baresys]
		mset2 = MembraneSet()
		c0trash = mset2.load_points_vtu(vtudirb,extra_props='induced_cur',start=startb,end=endb,nbase=nbase,
			lenscale=lenscale)
		mset2.surfacer()
		
#---calculate mode couplings
if int(seq[1]):
	mset.calculate_undulation_spectrum(lenscale=1.)
	mset.analyze_undulations(lenscale=1.)
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
	extraterms = -1*t1*qmagst**2-t2*qmagst**2+t3
	termsumspec = array([[qmagst[i,j],termsum[i,j]] for j in range(nt) for i in range(mt)])
	t0spec2 = array([[qmagst[i,j],uqrawmean[i,j]*qmagst[i,j]**4] for j in range(nt) for i in range(mt)])
	t1spec2 = array([[qmagst[i,j],t1[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t2spec2 = array([[qmagst[i,j],t2[i,j]*qmagst[i,j]**2] for j in range(nt) for i in range(mt)])
	t3spec2 = array([[qmagst[i,j],t3[i,j]] for j in range(nt) for i in range(mt)])
	if barecompare:
		mset2.calculate_undulation_spectrum(lenscale=1.)
		mset2.analyze_undulations(lenscale=1.)
		grid = mset2.griddims
		m,n = grid[0]-1,grid[1]-1
		hqs_bare = [fftwrap(mset2.surf[i])/double(m*n) for i in range(len(mset2.surf))]
		hqsa_bare = autocorr(hqs_bare,direct=0)
		t0_bare = uqrawmean = mean((abs(array(hqsa_bare)))**2,axis=0)
		cqs = [fftwrap(c0s[i])/double(m*n) for i in range(len(c0s))]
		Lx,Ly = mean(mset2.vecs,axis=0)[0:2]
		cm,cn = [int(round(i/2.-1)) for i in shape(hqs_bare[0])]
		qmags2 = array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) 
			for j in range(0,n)] for i in range(0,m)])
		qmagfilter=[10**-10,10**6]
		spectrum = array(sorted(zip(*[mset2.qrawmean,mset2.uqrawmean,mset2.uqrawstd]),
			key=lambda x: x[0]))[1:]
		spectrumf = array(filter(lambda x: x[0] >= qmagfilter[0] and x[0] <= qmagfilter[1], spectrum))
		specbare = array(filter(lambda x: x[0] >= 1./10.*qmagfilter[0] and 
			x[0] <= qmagfilter[1]*10.,spectrum))
		energy_bare = array([[specbare[j,0],specbare[j,0]**4*specbare[j,1]] for j in range(len(specbare))])
	
#---plots
if int(seq[2]):

	#---settings
	islognorm = True
	vmin = 10**-10
	vmax = 10**0
	cmap = mpl.cm.jet

	#---figure for 2D plots
	fig = plt.figure(figsize=(8,8))
	gs = gridspec.GridSpec(2,2)

	#---2D plots
	axes = []
	ax = plt.subplot(gs[0,0])
	axes.append(ax)
	ax.set_title('term 1')
	ax.imshow(array(t0).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$',fontsize=18)
	ax = plt.subplot(gs[0,1])
	axes.append(ax)
	ax.set_title('term 2')
	im1 =ax.imshow(array(t1).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}C_{0,-q}\right\rangle}$',fontsize=18)
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(im1,cax=cax1)
	ax = plt.subplot(gs[1,0])
	axes.append(ax)
	ax.set_title('term 3')
	ax.imshow(array(t2).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle C_{0,q}h_{-q}\right\rangle}$',fontsize=18)
	ax = plt.subplot(gs[1,1])
	axes.append(ax)
	ax.set_title('term 4')
	im3 = ax.imshow(array(t3).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		cmap=cmap,norm=(mpl.colors.LogNorm() if islognorm else None),vmin=vmin,vmax=vmax)
	ax.set_title(r'$\mathbf{\left\langle C_{0,q}C_{0,-q}\right\rangle}$',fontsize=18)
	divider = make_axes_locatable(ax)
	cax3 = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(im3,cax=cax3)

	for a in range(len(axes)):
		ax = axes[a]
		if a != 0:
			ax.set_yticks([])
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])	
		if a in [0,2]:
			ax.set_ylabel(r'$\left|\mathbf{q_y}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		if a in [2,3]:
			ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm^{-1}})$',fontsize=14)

	plt.savefig(pickles+'fig-bilayer-couple-view-terms-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()
	plt.cla()
	plt.clf()

#---plots, compare undulation spectra between bare and protein systems
if int(seq[3]) and barecompare:

	#---figure for 2D plots
	fig = plt.figure(figsize=(8,8))
	gs2 = gridspec.GridSpec(1,2,wspace=0.4)
	axes = []
	m,n = mset.griddims
	imdat0 = (abs(mean(mset.uqraw,axis=0))/double(m*n)**2)
	m,n = mset2.griddims
	imdat1 = (abs(mean(mset2.uqraw,axis=0))/double(m*n)**2)
	extrema = [min(imdat0.min(),imdat1.min()),max(imdat0.max(),imdat1.max())]

	ax = plt.subplot(gs2[0])
	axes.append(ax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$, protein',fontsize=18)
	im = ax.imshow(array(imdat0).T, extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,norm=mpl.colors.LogNorm(),vmin=extrema[0],vmax=extrema[1])

	ax = plt.subplot(gs2[1])
	axes.append(ax)
	im = ax.imshow(array(imdat1).T,extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,norm=mpl.colors.LogNorm(),vmin=extrema[0],vmax=extrema[1])
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$, bare',fontsize=18)
	if 0:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(im,cax=cax)

	for a in range(len(axes)):
		ax = axes[a]
		ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		if a in [0]:		
			ax.set_ylabel(r'$\left|\mathbf{q_y}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm+1)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn+1)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])

	plt.savefig(pickles+'fig-bilayer-couple-view-comparehqhq-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()
	plt.cla()
	plt.clf()

#---plots, compare undulation spectra between bare and protein systems scaled by q4
if int(seq[4]) and barecompare:

	#---figure for 2D plots
	fig = plt.figure(figsize=(8,8))
	gs2 = gridspec.GridSpec(1,2,wspace=0.4)
	axes = []
	m,n = mset.griddims
	imdat0 = t0*qmagst**4
	m,n = mset2.griddims
	imdat1 = t0_bare*qmagst**4
	extrema = [min(imdat0.min(),imdat1.min()),max(imdat0.max(),imdat1.max())]

	ax = plt.subplot(gs2[0])
	axes.append(ax)
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle {\left|\mathbf{q_y}\right|}^{4}}$, protein',fontsize=18)
	im = ax.imshow(array(imdat0).T, extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,vmin=extrema[0],vmax=extrema[1])

	ax = plt.subplot(gs2[1])
	axes.append(ax)
	im = ax.imshow(array(imdat1).T,extent=None,interpolation='nearest',
		aspect='equal',origin='lower',cmap=cmap,vmin=extrema[0],vmax=extrema[1])
	ax.set_title(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle {\left|\mathbf{q_y}\right|}^{4}}$, bare',fontsize=18)
	if 0:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(im,cax=cax)

	for a in range(len(axes)):
		ax = axes[a]
		ax.set_xlabel(r'$\left|\mathbf{q_x}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		if a in [0]:			
			ax.set_ylabel(r'$\left|\mathbf{q_y}\right|(\mathrm{nm^{-1}})$',fontsize=14)
		ax.set_xticks(array(list(arange(0,m/2,sskip)*-1)[:0:-1]+list(arange(0,m/2,sskip)))+cm)
		ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskip)*-1)[:0:-1]+
			list(arange(0,m/2,sskip)))*lenscale])
		ax.set_yticks(array(list(arange(0,n/2,sskip)*-1)[:0:-1]+list(arange(0,n/2,sskip)))+cn)
		ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskip)*-1)[:0:-1]+
			list(arange(0,n/2,sskip)))*lenscale])

	plt.savefig(pickles+'fig-bilayer-couple-view-comparehqhqq4-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()
	plt.cla()
	plt.clf()

#---plots
if int(seq[5]):

	#---separate 1D plot
	fig = plt.figure(figsize=(12,6))
	gs3 = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0)

	#---combined 1D spectra
	ax = plt.subplot(gs3[0])

	#---plot undulation spectra 1D, code cribbed from plotter.py
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.grid(True)
	spectrum = array(sorted(zip(*[mset.qrawmean,mset.uqrawmean,mset.uqrawstd]), key=lambda x: x[0]))[1:]
	ax.scatter(spectrum[:,0],spectrum[:,1],marker='o',color='k')
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
	ax.plot(xmod2,ymod2,color='k',linestyle='-',linewidth=1.5)
	ax.text(0.05,0.05,r'$\mathbf{\kappa = %3.2f}$'%kappa, transform=ax.transAxes,fontsize=14)
	ax.set_ylabel(r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$',fontsize=14)
	ax.scatter(spectrum1d[:,0],spectrum1d[:,1],color=clrs[0],
		label=r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$')
	ax.scatter(t1spec[:,0],t1spec[:,1],color=clrs[1],marker='+')	
	ax.scatter(t2spec[:,0],t2spec[:,1],color=clrs[2],marker='x',
		label=r'$\mathbf{\left\langle C_{0,q}h_{-q}\right\rangle}$')
	ax.scatter(t3spec[:,0],t3spec[:,1],color=clrs[3],marker='.',
		label=r'$\mathbf{\left\langle C_{0,q}C_{0,-q}\right\rangle}$')
		
	if barecompare:
		ax.scatter(specbare[:,0],specbare[:,1],color=clrs[4],marker='.',
			label=r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle}$, bare')
	ax.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm^{-1}})$',fontsize=14)
	ax.set_title('undulations')
	plt.legend(loc='upper right',fontsize=10)

	ax = plt.subplot(gs3[1:])
	reord0 = np.lexsort((termsumspec[:,1],termsumspec[:,0]))
	ax.plot(termsumspec[reord0,0],termsumspec[reord0,1],'.-',color='k',label='energy',lw=2)
	corrects = array([[termsumspec[i,0],t0spec2[i,1]] 
		for i in range(len(termsumspec))])
	reord1 = np.lexsort((corrects[:,1],corrects[:,0]))
	ax.plot(termsumspec[reord1,0],corrects[reord1,1],'.-',color=clrs[0],lw=2,
		label=r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle} {\left|\mathbf{q}\right|}^4$')
	if barecompare:
		reord2 = np.lexsort((energy_bare[:,1],energy_bare[:,0]))
		ax.plot(energy_bare[reord2,0],energy_bare[reord2,1],'.-',lw=1.5,color=clrs[4],
		label=r'$\mathbf{\left\langle h_{q}h_{-q}\right\rangle} {\left|\mathbf{q}\right|}^4$, bare')	
	if 0:
		ax.fill_between(termsumspec[reord1,0],termsumspec[reord1,1],
			energy_bare[reord2,1],facecolor='k',alpha=0.35)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.grid(True)
	ax.yaxis.tick_right()
	ax.set_ylabel(r'$\left\langle \mathscr{H}_{{\rm el}}\right\rangle \left[\kappa{\cal L}^{2}k_{B}T\right]^{-1}$',fontsize=18)
	ax.set_xlabel(r'$\left|\mathbf{q}\right|(\mathrm{nm^{-1}})$',fontsize=14)
	ax.yaxis.set_label_position("right")
	if 1:
		ax.set_ylim((0.4*min(termsumspec[:,1]),1.1*max(max(termsumspec[:,1]),
			max(corrects[:,1]),max(energy_bare[:,1]))))
	else:
		ax.set_ylim((10**-6,10**-4))
	ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
	plt.legend(loc='lower right')
	ax.set_title('energy')
	ax.grid(True, which='minor')
	plt.savefig(pickles+'fig-bilayer-couple-view-spectra1d-'+testname+'.png',
		dpi=500,bbox_inches='tight')
	if showplots:
		plt.show()
	plt.cla()
	plt.clf()
	

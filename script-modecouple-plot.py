#!/usr/bin/python

import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---settings
fspec = 'homog' if hypothesis['kappa']['fore'] == hypothesis['kappa']['back'] else 'heterog'
combosoln_rowsum = [
	range(-10,0),
	range(10),
	][1]
toscreen = False

#---MAIN
#-------------------------------------------------------------------------------------------------------------

def plot_spectrum(dat,denote,logplot=False,show=False):
	'''Temporary spectrum plotting program.'''
	datflat = abs(array(dat))+(10**-10 if logplot else 0)
	print datflat.max()
	print datflat.min()
	ax = plt.subplot(111)
	im = ax.imshow(
		datflat.T,
		interpolation='nearest',
		origin='lower',
		norm=(mpl.colors.LogNorm() if logplot else None))
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	if not show: 
		plt.savefig(
			allsets['pickles']+'FIGURES/fig-modecouple-'+denote+\
			('-log-1' if logplot else '')+'.png',dpi=300)
		plt.clf()
	else: plt.show()

#---combine some rows to look at the coupling
fullans = ms.fullans
n2,m2 = ms.m2,ms.n2

#---deprecated rowcombo
rowcombo = mean([array([[fullans[1][fr][j+(i-1)*m2] for j in range(n2)] 
	for i in range(m2)])*fullans[0][fr]**2 for fr in combosoln_rowsum],axis=0)
raws = [fullans[0][i]*reshape(fullans[1][i],(m2,n2)) for i in combosoln_rowsum]
combosoln = mean(raws,axis=0)

if 0:
	#---loop over plots
	for logplot in [True,False][:1]:
		for data in [
			[ms.kqqp,'kqqp'],
			[ms.kfield,'kfield'],
			[ms.full_matrix,'full_matrix'],
			[fullans[1],'fullans'],
			[combosoln,'combosoln'],
			[rowcombo,'rowcombo'],
			][:]:
			print 'data = '+data[1]+' logplot = '+str(logplot)
			plot_spectrum(data[0],callsign+'-term-'+data[1]+'-'+fspec,logplot=logplot,show=toscreen)

	#---plot eigenspectrum
	if 0:
		ax=plt.subplot(111);
		ax.plot(range(len(fullans[0]))[:100],abs(fullans[0])[:100]);
		ax.set_xscale('log');
		ax.set_yscale('log');plt.show()
if 1: 
	plot_spectrum(real((ms.qmags*ms.qmags)*(ms.qmags*ms.qmags)*termlist[0]).T,'',show=True,logplot=True)
	plot_spectrum(real(termlist[0]).T,'',show=True,logplot=True)
if 0:
	if 0: plot_spectrum(xvals,'',show=True,logplot=True)
	yvals = sum(abs(ms.full_matrix),axis=1)[1:]
	if 0: xvals = linalg.norm(ms.usshift,axis=1)[1:]
	xvals = ms.qmags[0][1:]
	ax = plt.subplot(111)
	ax.plot(xvals,yvals,'o')
	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim((-0.01,100))
	#ax.set_ylim((10**1,10**6))
	plt.show()
if 0:
	qmagshift = sqrt(sum((ms.usshift/(array([Lx,Ly])*ms.mset.lenscale/pi))**2,axis=1))


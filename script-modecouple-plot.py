#!/usr/bin/python

import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def tmpspecplotter(dat,denote,logplot=False,show=False):
	'''Temporary spectrum plotting program.'''
	datflat = abs(array(dat))+(10**-10 if logplot else 0)
	print datflat.max()
	print datflat.min()
	ax = plt.subplot(111)
	if logplot:
		im = ax.imshow(
			datflat.T,
			interpolation='nearest',
			origin='lower',
			norm=mpl.colors.LogNorm()
			)
	else:
		im = ax.imshow(
			datflat.T,
			interpolation='nearest',
			origin='lower',
			vmax=datflat.max(),
			vmin=datflat.min()
			)
	if 1:
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	if not show: plt.savefig(
		allsets['pickles']+'FIGURES/fig-modecouple-'+denote+\
		('-log-1' if logplot else '')+'.png',dpi=300)
	else: plt.show()
	plt.clf()

if 'fullans' not in globals() or 1:
	fullans = ms.fullans
	n2,m2 = ms.m2,ms.n2
	rowcombo = mean([array([[fullans[1][fr][j+(i-1)*m2] for j in range(n2)] 
		for i in range(m2)])*fullans[0][fr]**2 for fr in range(100)],axis=0)
if 1:
	for logplot in [False,True]:
		for data in [
			[ms.kqqp,'kqqp'],
			[ms.kfield,'kfield'],
			[ms.full_matrix,'full_matrix'],
			[fullans[1],'fullans'],
			[rowcombo,'rowcombo1'],
			][1:]:
			print 'data = '+data[1]+' logplot = '+str(logplot)
			tmpspecplotter(data[0],callsign+'-term-'+data[1]+'-homog',logplot=logplot,show=False)


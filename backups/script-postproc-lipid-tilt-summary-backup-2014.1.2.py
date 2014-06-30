#!/usr/bin/python -i

#---this one is for comparing distance to proteins

if 0:

	from membrainrunner import *

	location = ''
	execfile('locations.py')
	execfile('plotter.py')
	
	tilt_pickles = ['pkl.tilt.v700.md.part0002.100000-200000-200.pkl',
		'pkl.tilt.v701.md.part0003.60000-160000-200.pkl',
		'pkl.tilt.membrane-v550.md.part0006.200000-300000-200.pkl']

	#mset = unpickle(pickles+tilt_pickles[-1])
	
	allangles = []
	sysnames = []
	msets = []
	for i in range(len(tilt_pickles)):
		mset = unpickle(pickles+tilt_pickles[i])
		msets.append(mset)
		sysnames.append(tilt_pickles[-2][18:-4])
		allangles.append([180.*k for k in array(flatten(mset.store[0].get(['monolayer',0,'type','angle'])))])
		
	
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

'''
...
monitor
'''
#---plot the raw histogram
if 0:
	names = ['parallel','antiparallel','control']
	clrs = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
	hists = []
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.set_xlim((40,180))
	for i in range(len(angles)):
		ang = angles[i]
		hist0 = numpy.histogram(ang,bins=100,normed=True)
		hists.append(hist0)
		#plt.bar(hist0[1][:-1],2*hist0[0],color='b',alpha=0.5)
		c = clrs[i%len(clrs)]
		plt.plot(hist0[1][1:],hist0[0],'-',label=names[i],lw=2,color=c)
	#ax.set_yscale('log')
	plt.xlabel('Angle', labelpad = 10)
	plt.ylabel('Frequency', labelpad = 10)
	plt.title('lipid tilt angles')
	plt.legend(loc=2)
	plt.show()
	
#---plot the angle vs distance histograms
if 1:
	mset = msets[0]
	dists = mset.store[1].get(['type','mindist','monolayer',0])
	angles = mset.store[0].get(['monolayer',0,'type','angle'])
	distsflat = [i for j in dists for i in j]
	anglesflat = [180.*k for k in [i for j in angles for i in j]]
	mset = msets[1]
	dists = mset.store[1].get(['type','mindist','monolayer',0])
	angles = mset.store[0].get(['monolayer',0,'type','angle'])
	distsflat2 = [i for j in dists for i in j]
	anglesflat2 = [180.*k for k in [i for j in angles for i in j]]
	yspan = (min(min(anglesflat),min(anglesflat2)),max(max(anglesflat),max(anglesflat2)))
	xspan = (min(min(distsflat),min(distsflat2)),max(max(distsflat),max(distsflat2)))
	
	fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(14,7))
	ticknums = 6

	H, xedges, yedges = histogram2d(distsflat,anglesflat,range=(xspan,yspan),bins=50)
	midx = (xedges[1:]+xedges[:-1])/2
	midy = (yedges[1:]+yedges[:-1])/2
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	cmap.set_bad(cmap(0),1.)
	normH = [[float(H[i][j])/sum(H,axis=0)[j] for j in range(50)] for i in range(50)]
	axes[0][0].imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=None,cmap=cmap)
	xts = range(int(round(min(midx),-1)),int(round(max(midx),-1)),int(round((max(midx)-min(midx))/ticknums,-1)))
	axes[0][0].axes.set_xticks([[round(i,-1) for i in midx].index(j) for j in xts])
	axes[0][0].axes.set_xticklabels(xts)
	yts = range(int(round(min(midy),-1)),int(round(max(midy),-1)),int(round((max(midy)-min(midy))/ticknums,-1)))
	axes[0][0].axes.set_yticks([[round(i,-1) for i in midy].index(j) for j in yts])
	axes[0][0].axes.set_yticklabels(yts)
	#axes[0][0].axes.set_xlabel('minimum distance to protein (\AA)')
	axes[0][0].axes.set_ylabel('tilt angle (degrees)')
	axes[0][0].axes.set_title('parallel')

	H, xedges, yedges = histogram2d(distsflat2,anglesflat2,bins=50)
	midx = (xedges[1:]+xedges[:-1])/2
	midy = (yedges[1:]+yedges[:-1])/2
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	cmap.set_bad(cmap(0),1.)
	normH = [[(float(H[i][j])/sum(H,axis=0)[j] if H[i][j] != 0 else 10**-10) for j in range(50)] for i in range(50)]
	axes[0][1].imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=None,cmap=cmap)
	xts = range(int(round(min(midx),-1)),int(round(max(midx),-1)),int(round((max(midx)-min(midx))/ticknums,-1)))
	axes[0][1].axes.set_xticks([[round(i,-1) for i in midx].index(j) for j in xts])
	axes[0][1].axes.set_xticklabels(xts)
	yts = range(int(round(min(midy),-1)),int(round(max(midy),-1)),int(round((max(midy)-min(midy))/ticknums,-1)))
	axes[0][1].axes.set_yticks([[round(i,-1) for i in midy].index(j) for j in yts])
	axes[0][1].axes.set_yticklabels(yts)
	#axes[0][0].axes.set_xlabel('minimum distance to protein (\AA)')
	#axes[0][1].axes.set_ylabel('tilt angle (degrees)')
	axes[0][1].axes.set_title('antiparallel')


	H, xedges, yedges = histogram2d(distsflat,anglesflat,range=(xspan,yspan),bins=50)
	midx = (xedges[1:]+xedges[:-1])/2
	midy = (yedges[1:]+yedges[:-1])/2
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	cmap.set_bad(cmap(0),1.)
	normH = [[float(H[i][j])/sum(H,axis=0)[j] for j in range(50)] for i in range(50)]
	axes[1][0].imshow(array(normH).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',norm=None,cmap=cmap,vmin=0,vmax=0.1)
	xts = range(int(round(min(midx),-1)),int(round(max(midx),-1)),int(round((max(midx)-min(midx))/ticknums,-1)))
	axes[1][0].axes.set_xticks([[round(i,-1) for i in midx].index(j) for j in xts])
	axes[1][0].axes.set_xticklabels(xts)
	yts = range(int(round(min(midy),-1)),int(round(max(midy),-1)),int(round((max(midy)-min(midy))/ticknums,-1)))
	axes[1][0].axes.set_yticks([[round(i,-1) for i in midy].index(j) for j in yts])
	axes[1][0].axes.set_yticklabels(yts)
	axes[1][0].axes.set_xlabel('minimum distance to protein (\AA)')
	axes[1][0].axes.set_ylabel('tilt angle (degrees)')
	axes[1][0].axes.set_title('parallel (normed)')

	H, xedges, yedges = histogram2d(distsflat2,anglesflat2,bins=50)
	midx = (xedges[1:]+xedges[:-1])/2
	midy = (yedges[1:]+yedges[:-1])/2
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	cmap.set_bad(cmap(0),1.)
	normH = [[float(H[i][j])/sum(H,axis=0)[j] for j in range(50)] for i in range(50)]
	axes[1][1].imshow(array(normH).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=None,cmap=cmap,vmin=0,vmax=0.1)
	xts = range(int(round(min(midx),-1)),int(round(max(midx),-1)),int(round((max(midx)-min(midx))/ticknums,-1)))
	axes[1][1].axes.set_xticks([[round(i,-1) for i in midx].index(j) for j in xts])
	axes[1][1].axes.set_xticklabels(xts)
	yts = range(int(round(min(midy),-1)),int(round(max(midy),-1)),int(round((max(midy)-min(midy))/ticknums,-1)))
	axes[1][1].axes.set_yticks([[round(i,-1) for i in midy].index(j) for j in yts])
	axes[1][1].axes.set_yticklabels(yts)
	axes[1][1].axes.set_xlabel('minimum distance to protein (\AA)')
	#axes[1][1].axes.set_ylabel('tilt angle (degrees)')
	axes[1][1].axes.set_title('antiparallel (normed)')

	H, xedges, yedges = histogram2d(distsflat2,anglesflat2,bins=50)
	midx = (xedges[1:]+xedges[:-1])/2
	midy = (yedges[1:]+yedges[:-1])/2
	extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
	cmap = mpl.cm.jet
	cmap.set_bad(cmap(0),1.)
	normH = [[float(H[i][j])/sum(H,axis=0)[j] for j in range(50)] for i in range(50)]
	hist0 = numpy.histogram([180*k for k in angles[-1]],bins=50,normed=True,range=yspan)
	axes[0][2].imshow(array([list(hist0[0]) for i in range(50)]).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=None,cmap=cmap,vmin=0,vmax=0.1)
	xts = range(int(round(min(midx),-1)),int(round(max(midx),-1)),int(round((max(midx)-min(midx))/ticknums,-1)))
	axes[0][2].axes.set_xticks([[round(i,-1) for i in midx].index(j) for j in xts])
	axes[0][2].axes.set_xticklabels(xts)
	yts = range(int(round(min(midy),-1)),int(round(max(midy),-1)),int(round((max(midy)-min(midy))/ticknums,-1)))
	axes[0][2].axes.set_yticks([[round(i,-1) for i in midy].index(j) for j in yts])
	axes[0][2].axes.set_yticklabels(yts)
	axes[0][2].axes.set_title('control')

	names = ['parallel','antiparallel','control']
	clrs = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
	hists = []
	axes[1][2].set_xlim((40,180))
	for i in range(len(allangles)):
		ang = allangles[i]
		hist0 = numpy.histogram(ang,bins=100,normed=True)
		hists.append(hist0)
		c = clrs[i%len(clrs)]
		axes[1][2].plot(hist0[1][1:],hist0[0],'-',label=names[i],lw=2,color=c)
	plt.xlabel('Angle', labelpad = 10)
	plt.ylabel('Frequency', labelpad = 10)
	axes[1][2].legend(loc=2)	
	
	fig.tight_layout()
	plt.savefig(pickles+'fig-lipid-tilts-summary-v700.v701.v550'+'.png', dpi=300)
	plt.show()
	


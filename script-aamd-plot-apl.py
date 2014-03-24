#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

if 'msets' not in globals(): 
	msets = []
	msets.append(unpickle(pickles+\
		'pkl.lipidarea.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl'))
mset = msets[0]
dat = mset.store[0]

clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in range(9)]
clrset = [clrs[i] for i in [8,7,1,0]]

def cellplot(dat,fig,altdat=None,vmin=None,vmax=None,fr=None,panels=1):
	rows,cols = shape(dat)
	gs = gridspec.GridSpec(rows,cols,wspace=0.0,hspace=0.0)
	#---rows are monolayers, columns are systems
	for row in range(rows):
		for col in range(cols):
			vor = dat[row][col]
			ax = fig.add_subplot(gs[row,col],aspect='equal')
			regionlist = [vor.regions[i] for i in vor.point_region]
			for r in range(len(regionlist)):
				region = regionlist[r]
				if not -1 in region:
					polygon = [vor.vertices[i] for i in region]
					if 0: axes.fill(*zip(*polygon),alpha=0.5)
					p = mpl.patches.Polygon(polygon,alpha=0.65,
						facecolor=('w' if r >= len(colorcodes[row]) else clrset[colorcodes[row][r]]),
						lw=0.5,edgecolor='k')
					ax.add_patch(p)
			for p in range(len(vor.points)):
				pt = vor.points[p]
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1.,facecolor=rand_color_list[p],
					edgecolor='k',lw=0.3))
			ax.set_xlim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][0])
			ax.set_ylim([[-0.1*i,1.1*i] for i in mean(mset.vecs,axis=0)[:2]][1])
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xticks([])
			ax.set_yticks([])
	return [gs,ax]

#---relative resids
resids = [[mset.resids_reracker.index(i) for i in j] for j in mset.resids]
#---color codes
colorcodes = [[[i for i in range(4) if j in resids[i]][0] for j in mono] for mono in mset.monolayer_residues]

#---populate the data structure for the movie maker
dat_vor = [[0. for cols in range(1)] for rows in range(2)]
#---rows are monolayers, columns are systems
dat_vor[0][0] = list(dat.get(['monolayer',0,'type','voronoi']))
dat_vor[1][0] = list(dat.get(['monolayer',1,'type','voronoi']))
rand_color_list = [np.random.rand(3,1) for i in range(2*len(dat_vor[0][0][0].points))]

plotmov(dat_vor,'celltest',plotfunc='cellplot')




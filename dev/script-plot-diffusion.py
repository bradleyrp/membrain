#!/usr/bin/python -i

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

from membrainrunner import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

#---Analysis parameters
skip = None
framecount = None
location = ''
execfile('locations.py')

#---selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---analysis plan
analysis_descriptors = [
	(['membrane-v530'],'resname PI2P and name P',director_asymmetric)]
analyses = [analysis_descriptors[i] for i in range(len(analysis_descriptors))]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def make_colormap(seq):
	'''
	Cribbed from stack overflow.
	Return a LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	'''
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = seq[i - 1]
			r2, g2, b2 = seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)

#---MAIN
#-------------------------------------------------------------------------------------------------------------

ad = analyses[-1]
(tests,selector,director) = ad
mset = MembraneSet()

#---custom
basedir = '/home/rpb/tmp/'
gro = 'system-input.pi2p.gro'
traj = 'system.pi2p.xtc'

#---prepare
mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
mset.identify_monolayers(director)
allselect_lipids = mset.universe.selectAtoms(selector)
for mononum in range(2):
	validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
		set([i-1 for i in allselect_lipids.resids()])))
	mset.selections.append(sum([allselect_lipids.residues[
		list(allselect_lipids.resids()).index(i+1)].selectAtoms(selector) 
		for i in validresids]))

#---frame selection header
end = None
start = None
if framecount == None:
	if end == None: end = mset.nframes
	if start == None: start = 0
	if skip == None: skip = 1
else:
	start = 0
	end = mset.nframes
	skip = int(float(mset.nframes)/framecount)
	skip = 1 if skip < 1 else skip
	
#---collect points
lipidpts = []
for frameno in range(start,end,skip):
	vecs = mset.vec(frameno)
	pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
	lipidpts.append(pts1)
	
#---plot
c = mcolors.ColorConverter().to_rgb
vecs=mean(mset.vecs,axis=0)
plotslice = slice(None,None,1)
selected_lipids = [i for i in range(0,40)]
fig = plt.figure(figsize=(8,8))	
ax = fig.add_subplot(111,aspect='equal')
ax = plt.subplot(111)
ax.set_title('diffusion map')
prev_rando = 0
for sel in selected_lipids:
	print sel
	rando = random.randint(0,len(selected_lipids)) + prev_rando
	cmap = make_colormap([c('black'),0.5,mpl.cm.jet(float(int(rando)%len(selected_lipids))/\
		len(selected_lipids))[:-1]])
	alphas = np.abs(np.linspace(-1.0, 1.0,cmap.N))
	cmap._init()
	cmap._lut[:-3,-1] = alphas
	col = [cmap(float(i)/(len(lipidpts[0]))) for i in xrange(len(lipidpts[0]))]
	#ax.scatter(array(lipidpts)[plotslice,sel,0],array(lipidpts)[plotslice,sel,1],
	#	marker='.',lw = 2,s=20,c=col,edgecolors='none')
	ax.plot(array(lipidpts)[plotslice,sel,0],array(lipidpts)[plotslice,sel,1],
		marker='.',lw = 1)
ax.set_xlim((-0.2*vecs[0],1.2*vecs[0]))
ax.set_ylim((-0.2*vecs[1],1.2*vecs[1]))
ax.axvline(x=0,linewidth=2, color='k')
ax.axvline(x=vecs[0],linewidth=2, color='k')
ax.axhline(y=0,linewidth=2, color='k')
ax.axhline(y=vecs[1],linewidth=2, color='k')
ax.set_xticks(arange(-0.2*vecs[0],1.2*vecs[0],10))
ax.set_yticks(arange(-0.2*vecs[1],1.2*vecs[1],10))
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.grid(True, which='minor')
plt.show()


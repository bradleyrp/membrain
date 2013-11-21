#!/bin/python

#import curvcalc
#import numpy

self = mset
points = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
	for i in self.monolayer_residues[0]])
tmp = mset.tri[0]
protcom = mean(mset.protein[0],axis=0)
subji = numpy.argmin([norm(i[0:2]-protcom[0:2]) for i in points])
#subji = tmp.simplices[int(tmp.find_simplex(protcom[0:2]))][1]
#neighborsi = list(set([v for w in [list(tmp.simplices[i]) for i in tmp.neighbors[subji]] for v in w]))
neighborsi = list(set([w for v in [list(i) for i in tmp.simplices if subji in i] for w in v]))
nneighbors = len(neighborsi)
neighbors = array([points[i] for i in neighborsi])
selfpts = points[subji]
print selfpts
print nneighbors
print neighbors
outpts = numpy.concatenate(([selfpts],neighbors),axis=0)

meshplot(points,show='wire')
meshplot(outpts)
meshpoints(protcom,scale_factor=20,color=(1,1,1))
meshpoints(selfpts,scale_factor=20,color=(1,1,1))
#meshpoints(mset.protein[0],scale_factor=20,color=(1,1,1))

#curvcalc(len(outpts),outpts)


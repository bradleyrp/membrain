#!/usr/bin/python -i

from membrainrunner import *
import numpy
from scipy.spatial import Delaunay
import memint
location = ''
execfile('locations.py')
execfile('plotter.py')

mset = unpickle(pickles+'pkl.structures.membrane-v700.md.part0006.360000-460000-200.pkl')
points = mset.unzipgrid(mset.surf[0],vecs=mset.vec(0))
dl = Delaunay(points[:,0:2])
#print '\ncalculating neighborlist'
#nlist = [where([any(i==j) for i in dl.simplices])[0] for j in range(len(points))]
#nlist = [[s for s in range(len(dl.simplices)) if any(dl.simplices[s]==i)] for i in range(len(points))]
#print 'done'
#protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]
#select = [i for i in range(len(points)) if norm(points[i]-protcom) < 400]
#sel_tris = [dl.simplices[i] for i in select]
memint.init(len(points),len(dl.simplices))
for p in range(len(points)):
	memint.addvertex(p,points[p],nlist[p])
for t in range(len(dl.simplices)):
	memint.addtriangle(t,dl.simplices[t])

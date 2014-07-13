#!/usr/bin/python -i

from membrainrunner import *
import numpy
#sys.path.append('/home/rpb/libs/scipy-0.13/lib64/python2.7/site-packages/scipy')
from scipy.spatial import Delaunay
#import imp
#spatial2 = imp.load_source('spatial', '/home/rpb/libs/scipy-0.13/lib64/python2.7/site-packages/scipy13/spatial/__init__.py')

import memint
execfile('plotter.py')

#---Analysis parameters
location = 'dark'
skip = 1
framecount = None

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-repository/'
	locations = '/home/rpb/worker/membrain/trajectory-rpb-light'	
	pickles = '/home/rpb/worker-big/membrane-repository/pickle-repository/'
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/store-delta/worker/worker-big/membrane-repository/'
	locations = '/store-delta/worker/membrain/trajectory-rpb-dark'
	pickles = '/store-delta/worker/worker-big/membrane-repository/pickle-repository/'
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

#...

mset = unpickle(pickles+'pkl.avgstruct.membrane-v612.md.part0003.skip10.pkl')
points = mset.unzipgrid(mset.surf[0],vecs=mset.vec(0))
dl = Delaunay(points[:,0:2])

(indices, indptr) = dl.vertex_neighbor_vertices
nlist = [list(indptr[indices[k]:indices[k+1]]) for k in range(dl.npoints)]
assoctris = [list(where(all([list(any([list(any(dl.simplices==j,axis=1)) for j in nlist[t0]],axis=0)),list(any(dl.simplices==t0,axis=1))],axis=0))[0]) for t0 in range(dl.npoints)]

protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]
#select = [i for i in range(len(points)) if norm(points[i]-protcom) < 400]
#sel_tris = [dl.simplices[i] for i in select]
memint.init(len(points),len(dl.simplices))
for p in range(len(points)):
	memint.addvertex(p,list(points[p]),nlist[p],assoctris[p])
for t in range(len(dl.simplices)):
	memint.addtriangle(t,array(dl.simplices[t]))
'''	
for p in range(len(points)):
	memint.getprop(part='vertex',subpart='coords',index=p)	
for t in range(len(dl.simplices)):
	memint.getprop(part='triangle',subpart='vert',index=t)
'''	
memint.check()
'''
for p in range(len(points)):
	memint.getprop(part='vertex',subpart='mcur',index=p)	
	memint.getprop(part='vertex',subpart='pdir',index=p)	
'''

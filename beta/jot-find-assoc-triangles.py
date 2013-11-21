#!/usr/bin/bash

#This is my attempt to find triangles associated with a particular vertex in Delaunay faster

'''

trying to get associated triangles faster
list(where(any([list(any(dl.simplices==j,axis=1)) for j in nlist[40]],axis=0))[0])

list(any([list(any(dl.simplices==j,axis=1)) for j in nlist[40]],axis=0))
list(any(dl.simplices==40,axis=1))

list(where(all([list(any([list(any(dl.simplices==j,axis=1)) for j in nlist[40]],axis=0)),list(any(dl.simplices==40,axis=1))],axis=0))[0])

fast?
assoctris = [list(where(all([list(any([list(any(dl.simplices==j,axis=1)) for j in nlist[t0]],axis=0)),list(any(dl.simplices==t0,axis=1))],axis=0))[0]) for t0 in range(dl.npoints)]

slow,previous?
assoctris = [[s for s in range(len(dl.simplices)) if any(dl.simplices[s]==j)] for j in range(len(points))]

results in minutes:

calculating associated triangles
done
0.314769363403

calculating associated triangles
done
0.0246350328128


previous neighborlists calculated by

#nlist = [where([any(i==j) for i in dl.simplices])[0] for j in range(len(points))]
#nlist = [[s for s in range(len(dl.simplices)) if any(dl.simplices[s]==i)] for i in range(len(points))]


but now scipy 0.13 has a faster method

'''

#!/usr/bin/python

import scipy.interpolate

#---Check the bin offsets
if 0:
	meshplot([[int(avginterpsurf[i,j]) for i in range(shape(avginterpsurf)[1])] for j in range(shape(avginterpsurf)[1])])

if 1:
	subj = mset.unzipgrid(array(results),vecs=vecs,grid=shape(results))
	grid=shape(results)
	subj = mset.unzipgrid(array(results),vecs=vecs,grid=grid)
	rbfobj = scipy.interpolate.Rbf(subj[:,0],subj[:,1],subj[:,2],epsilon=10.0,function='Linear',smooth=1000)
	ti1 = linspace(0,vecs[0],grid[0]);ti2 = linspace(0,vecs[1],grid[1])
	XI, YI = meshgrid(ti1, ti2)
	ZI = rbfobj(XI,YI)
	result2 = mset.unzipgrid(transpose(ZI)/10**2,vecs=vecs,reverse=0)
	meshplot(array(result2),vecs=vecs,show='surf')
	
if 0:
	grid=shape(results)
	resultsunzip = mset.unzipgrid(results/10**4,vecs=vecs)
	subj = mset.wrappbc(resultsunzip,vecs,mode='grow',growsize=0.2)
	#subj = mset.unzipgrid(array(results),vecs=vecs,grid=shape(results))
	rbfobj = scipy.interpolate.Rbf(subj[:,0],subj[:,1],subj[:,2],epsilon=1.0,function='Linear',smooth=1000)
	ti1 = linspace(0,vecs[0],grid[0]);ti2 = linspace(0,vecs[1],grid[1])
	XI, YI = meshgrid(ti1, ti2)
	ZI = rbfobj(XI,YI)
	result2 = mset.unzipgrid(transpose(ZI)/10**2,vecs=vecs,reverse=0)
	checkmesh(array(result2),vecs=vecs,show='surf',tess=4)


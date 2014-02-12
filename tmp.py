#!/bin/bash

self = mset
surfdata = c0sraw
skip=1
interp='best'
interpdata = []
fr = 0

print 'Interpolating splines, '+str(fr)
vecs = mset.vec(0)
grid = mset.griddims
newvecs = [vecs[i]-vecs[i]/(grid[i]-1) for i in range(2)]
newgrids = self.griddims[0]-1,self.griddims[1]-1

res0 = self.wrappbc(array([[self.xyzs[fr][j][0],self.xyzs[fr][j][1],surfdata[fr][j]] 
	for j in range(len(surfdata[fr]))]),vecs=newvecs,mode='grow')
vecs = self.vecs[fr]
steps = [vecs[i]/(grid[i]-1) for i in range(2)]
res1 = self.makemesh(res0,newvecs,newgrids,method=interp)
#rezip = self.rezipgrid(res1,diff=1)

#meshplot(res0)
#meshplot(res1)
#checkmesh(res1,vecs=newvecs,tess=4,show='surf')

xyz = res1
vecs=None
frameno=0
grid=None
rounder_vecs=[1.0,1.0]
reverse=0
diff=1
whichind=2

#original rezip
if 0:
	if grid == None: grid = self.griddims
	if vecs == None: vecs = self.vec(frameno)
	steps = [vecs[i]/(grid[i]-1) for i in range(2)] if diff == 0 else [vecs[i]/(grid[i]-1) for i in range(2)]
	poslookup = [[xyz[i][j]/steps[j] for j in range(2)] for i in range(len(xyz))]
	surfgrid = [[0. for i in range(grid[1]-1*(diff==1))] for j in range(grid[0]-1*(diff==1))]
	for i in range(len(xyz)): 
		#---Note: added the round command below changes calculate_midplane time from 0.25 to 0.35!
		if int(poslookup[i][0]) < grid[0] and int(poslookup[i][1]) < grid[1]:
			surfgrid[int(round(poslookup[i][0]))][int(round(poslookup[i][1]))] = xyz[i][whichind]
	rezip = surfgrid
	interpdata = [rezip]
#---new rezip
else:
	if grid == None and diff != 1: grid = self.griddims
	elif grid == None and diff == 1: grid = self.griddims[0]-1,self.griddims[1]-1
	if vecs == None: vecs = self.vec(frameno)
	print grid
	steps = [vecs[i]/(grid[i]-1) for i in range(2)] if diff == 0 else [vecs[i]/(grid[i]) for i in range(2)]
	print vecs
	print steps
	poslookup = [[xyz[i][j]/steps[j] for j in range(2)] for i in range(len(xyz))]
	surfgrid = [[0. for i in range(grid[1])] for j in range(grid[0])]
	for i in range(len(xyz)): 
		#---Note: added the round command below changes calculate_midplane time from 0.25 to 0.35!
		if int(poslookup[i][0]) < grid[0] and int(poslookup[i][1]) < grid[1]:
			surfgrid[int(round(poslookup[i][0]))][int(round(poslookup[i][1]))] = xyz[i][whichind]
	rezip = surfgrid
	interpdata = [rezip]

#!/usr/local/python

def rezipgrid2(self,xyz,vecs=None,frameno=0,grid=None,rounder_vecs=[1.0,1.0],reverse=0,diff=0,whichind=2):
	'''Turns a regular set of points in 3-space into a 2D matrix.'''
	if grid == None: grid = self.griddims
	if vecs == None: vecs = self.vec(frameno)
	steps = [vecs[i]/(grid[i]-1) for i in range(2)] if diff == 0 else [vecs[i]/(grid[i]-1) for i in range(2)]
	poslookup = [[xyz[i][j]/steps[j] for j in range(2)] for i in range(len(xyz))] 
	surfgrid = [[0. for i in range(grid[1]-1*(diff==1))] for j in range(grid[0]-1*(diff==1))]
	for i in range(len(xyz)): 
		#if int(poslookup[i][0]) < grid[0] and int(poslookup[i][1]) < grid[1]:
		surfgrid[int(poslookup[i][0])][int(poslookup[i][1])] = xyz[i][whichind]
	return surfgrid

self=mset
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
frameno=45
#do 45 if you are on light, 11 on dark
pbcadjust=1
rounder=4.0
interp='best'
storexyz=True

'''Find the midplane of a molecular dynamics bilayer.'''
lenscale = self.lenscale
self.gotoframe(frameno)
topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) for i in self.monolayer_residues[0]])
botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) for i in self.monolayer_residues[1]])
if storexyz:
	self.xyzs.append([topxyz,botxyz])
#---First we wrap PBCs
topxyzwrap = self.wrappbc(topxyz,vecs=self.vec(frameno),mode='grow')
botxyzwrap = self.wrappbc(botxyz,vecs=self.vec(frameno),mode='grow')
#---Triangulate the surface on a regular grid
topmesh = self.makemesh(topxyzwrap,self.vec(frameno),self.griddims,method=interp)
botmesh = self.makemesh(botxyzwrap,self.vec(frameno),self.griddims,method=interp)
#---Convert from points in 3-space to heights on a grid
self.monolayer1.append(topmesh)
self.monolayer2.append(botmesh)
#---Take the average surface
surfz = [[1./2*(rezipgrid2(self,topmesh,frameno=frameno)[i][j]+\
	rezipgrid2(self,botmesh,frameno=frameno)[i][j]) for i in range(self.griddims[0])] \
	for j in range(self.griddims[1])]
self.surf_position.append(mean(surfz))
surfz = surfz - mean(surfz)

whichind=2
vecs = None
grid = None
if grid == None: grid = self.griddims
if vecs == None: vecs = self.vec(frameno)
xyz=topmesh
steps = [vecs[i]/(grid[i]-1) for i in range(2)] if diff == 0 else [vecs[i]/(grid[i]-1) for i in range(2)]
#steps = [4.11029837,4.11029837]
print steps
poslookup = [[xyz[i][j]/steps[j] for j in range(2)] for i in range(len(xyz))] 
surfgrid = [[0. for i in range(grid[1]-1*(diff==1))] for j in range(grid[0]-1*(diff==1))]
for i in range(len(xyz)): 
	#if int(poslookup[i][0]) < grid[0] and int(poslookup[i][1]) < grid[1]:
	surfgrid[int(round(poslookup[i][0]))][int(round(poslookup[i][1]))] = xyz[i][whichind]
	
	
	
	
	
	
	

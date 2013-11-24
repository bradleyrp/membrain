#!/usr/local/python

self=mset
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
frameno=11
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
surfz = [[1./2*(self.rezipgrid(topmesh,frameno=frameno)[i][j]+\
	self.rezipgrid(botmesh,frameno=frameno)[i][j]) for i in range(self.griddims[0])] \
	for j in range(self.griddims[1])]
self.surf_position.append(mean(surfz))
surfz = surfz - mean(surfz)


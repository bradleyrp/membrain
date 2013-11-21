#!/bin/python

'''aborted attempt to plot the voronoi tesselation in 3d, which is tricky and waste of time'''

self = mset

if 0:
	mset.identify_monolayers(director,startframeno=3)
	mset.identify_residues(residues)
	#mset.triangulator('name P',framecount=10)

if 0:
	coords = mset.universe.selectAtoms('name CL').coordinates()
	coordst = array([mset.universe.residues[i].selectAtoms('name P').coordinates()[0] for i in mset.monolayer_residues[0]])
	coordsb = array([mset.universe.residues[i].selectAtoms('name P').coordinates()[0] for i in mset.monolayer_residues[1]])
	#meshpoints(coords,scale_factor=[5 for i in range(len(coords))],color=(1,1,1))
	#meshplot(coordst)
	#meshplot(coordsb)
	drawbox3d(mset.vecs[0])
	tmp = mset.store[1].data[0][0][0]

if 0:
	self = mset
	frameno = 1
	frame = self.universe.trajectory[frameno]
	if type(selector) == str:
		selector = [selector for i in range(len(self.resnames))]
	top = []
	for t in range(3):
		top.extend([mean(self.universe.residues[i].selectAtoms(selector[t]).coordinates(),axis=0) for i in self.monolayer_by_resid[0][t]])
	bot = []
	for t in range(3):
		bot.extend([mean(self.universe.residues[i].selectAtoms(selector[t]).coordinates(),axis=0) for i in self.monolayer_by_resid[1][t]])
	points = top
	points_pbc = mset.wrappbc(points,vecs=self.vec(frameno),mode='grow')
	vmap = scipy.spatial.Voronoi(points_pbc[:,0:3])
	#triangles = [[i for i in vmap.regions[vmap.point_region[p]]] for p in range(len(vmap.points))]
	#meshplot(vmap.points)
	points = bot
	points_pbc = mset.wrappbc(points,vecs=self.vec(frameno),mode='grow')
	vmap = scipy.spatial.Voronoi(points_pbc[:,0:3])
	#triangles = [[i for i in vmap.regions[vmap.point_region[p]]] for p in range(len(vmap.points))]
	#meshplot(vmap.points)

if 1:
	if type(selector) == str:
		selector = [selector for i in range(len(self.resnames))]
	top = []
	for t in range(3):
		top.extend([mean(self.universe.residues[i].selectAtoms(selector[t]).coordinates(),axis=0) for i in self.monolayer_by_resid[0][t]])
	frameno = mset.universe.trajectory.frame
	points = array(top)
	points_pbc = mset.wrappbc(points,vecs=self.vec(frameno),mode='grow')
	vmap = scipy.spatial.Voronoi(points_pbc[:,0:3])
	
	meshpoints([vmap.vertices[i] for i in allpts if i != -1])
	
	regs = [vmap.regions[i] for i in vmap.point_region]
	allpts = list(set([i for j in regs for i in j]))
	
	
	meshpoints(points,scale_factor=[5 for i in range(len(points))],color=(0,0,0))

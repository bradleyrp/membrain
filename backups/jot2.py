#!/bin/python

self=mset
selector = ['resname PI2P and name P','name NA']

if 1:
	mset.identify_monolayers(director,startframeno=3)
	mset.identify_residues(residues)
	#mset.triangulator('name P',framecount=10)
	#---Select atoms
	allselect_lipids = self.universe.selectAtoms(selector[0])
	allselect_ions = self.universe.selectAtoms(selector[1])
	validresids = list(set.intersection(set(self.monolayer_residues[0]),\
		set([i-1 for i in allselect_lipids.resids()])))
	sel1 = sum([allselect_lipids.residues[allselect_lipids.resids().index(i+1)].selectAtoms(selector[0]) \
		for i in validresids])
	validresids = list(set.intersection(set(self.monolayer_residues[1]),\
		set([i-1 for i in allselect_lipids.resids()])))
	sel2 = sum([allselect_lipids.residues[allselect_lipids.resids().index(i+1)].selectAtoms(selector[0]) \
		for i in validresids])
	self.selections.append(sel1)
	self.selections.append(sel2)
	self.selections.append(allselect_ions)

mset.gotoframe(1)

coords = mset.universe.selectAtoms('name NA').coordinates()
coordst = array([mset.universe.residues[i].selectAtoms('name P').coordinates()[0] for i in mset.monolayer_residues[0]])
coordsb = array([mset.universe.residues[i].selectAtoms('name P').coordinates()[0] for i in mset.monolayer_residues[1]])
meshplot(coordst)
meshplot(coordsb)
drawbox3d(mset.vecs[0])

#meshpoints(coords,scale_factor=[5 for i in range(len(coords))],color=(1,1,1))

frameno = 1
lipids_in = array(self.get_points(frameno,selection_index=0))
ions_in_all = self.get_points(frameno,selection_index=2)
bilayer_z = self.locate_bilayer(frameno)

#ions_in = array([c for c in ions_in_all if ((c[2]+-1*(int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*self.vec(frameno)[2])-bilayer_z) > 0.0])
#ions_in = array([c for c in ions_in_all if ((c[2]+-1*(int(c[2]/(self.vec(frameno)[2]/2+bilayer_z))%2)*self.vec(frameno)[2])-bilayer_z) > 0.0])
ions_in = array([c+[0,0,-1*(((c[2]>bilayer_z)*self.vec(frameno)[2])-0*bilayer_z)] for c in ions_in_all])

meshpoints(ions_in,scale_factor=[5 for i in range(len(ions_in))],color=(1,1,1))
meshpoints(ions_in_all,scale_factor=[3 for i in range(len(ions_in_all))],color=(0,0,0))

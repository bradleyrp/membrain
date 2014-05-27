#!/usr/bin/python

#---rewrap PBCs
self = mset
frameno = 10
mset.gotoframe(frameno)
topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) for i in self.monolayer_residues[0]])
topxyz2 = self.wrappbc(topxyz,vecs=self.vec(frameno),mode='grow')

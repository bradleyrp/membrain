#!/usr/bin/python

#---Fix the structure once, if necessary (ignore this rpb hack)
if 0:
	msets = [v509a,v509b]
	msets = [v514]
	for mset in msets:
		for i in range(len(mset.store)):
			mset.store[i].struct_opts = {'type' : {'points':0,'voronoi':1,'areas':2}}
			pickledump(mset,'pkl.'+mset.picklename+'.pkl')


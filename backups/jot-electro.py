#!/usr/bin/python
from MDAnalysis import *

if 0:
	fr = 10
	frame = mset.universe.trajectory[fr]
	pts = mset.universe.selectAtoms('all').coordinates()
	
if 0:
	tprpath = '/home/rpb/worker-big/membrane-repository/look-v510/md.part0033.tpr'
	structpath='/home/rpb/worker-big/membrane-repository/look-v510/system-input.gro'
	xtcpath='/home/rpb/worker-big/membrane-repository/look-v510//md.part0033.skip10.pbcmol.xtc'
	mset_get_charges = Universe(structpath,tprpath)

# import the pickle with charges
pickles = '/home/rpb/worker-big/membrane-repository/pickle-repository/'
mset_charges = unpickle(pickles+'pkl.v510.charges.pkl')


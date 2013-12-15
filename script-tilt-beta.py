#!/usr/bin/python -i

from membrainrunner import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = 1
framecount = 5
location = 'light'
execfile('locations.py')

#---Parameters
rounder = 10.0

#---Selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	(['membrane-v614'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None))]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

ad = analysis_descriptors[analysis_plan][0]
(tests,director,selector,protein_select,trajno) = ad
testno = 0
t = testno
traj = trajectories[systems.index(tests[t])][trajno][0]

'''...'''
mset = MembraneSet()
#---Load the trajectory
gro = structures[systems.index(tests[testno])]
basename = traj.split('/')[-1][:-4]
sel_surfacer = sel_cgmd_surfacer
print 'Accessing '+basename+'.'
mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
	resolution='cgmd')
#---Average structure calculation
mset.identify_monolayers(director,startframeno=0)
#mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,protein_selection=protein_select,end=150)

selector = selector_cgmd
self = mset

topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
	for i in self.monolayer_residues[0]])
botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
	for i in self.monolayer_residues[1]])


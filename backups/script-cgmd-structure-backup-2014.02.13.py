#!/usr/bin/python -i

from membrainrunner import *

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = None
location = ''
execfile('locations.py')

#---parameters
rounder = 20.0

#---selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
cgmd_protein = 'name BB'

#---analysis plan
analysis_plan = slice(None,None)
analysis_descriptors = [
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[360000,460000,200]),
	(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[60000,160000,200]),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[100000,200000,200]),
	(['membrane-v550'],director_cgmd,selector_cgmd,None,slice(-1,None),[300000,400000,200]),
	(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None),[500000,700000,400])]
analyses = [analysis_descriptors[i] for i in [0,-1,-2]]
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time()
print 'Starting analysis job.'
for ad in analyses:
	#---Load global variables with calculation specifications used in analysis functions above.
	(tests,director,selector,protein_select,trajno,timeslice) = ad
	for testno in range(len(tests)):
		print 'Running calculation: bilayer structure and undulations '+tests[testno]+'.'
		for traj in trajectories[systems.index(tests[testno])][trajno]:
			mset = MembraneSet()
			#---Load the trajectory
			gro = structures[systems.index(tests[testno])]
			basename = traj.split('/')[-1][:-4]
			sel_surfacer = sel_cgmd_surfacer
			print 'Accessing '+basename+'.'
			starttime = time.time()
			mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
				resolution='cgmd')
			print 'Load complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'
			#---Average structure calculation
			mset.identify_monolayers(director,startframeno=0)
			if protein_select == None:
				mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,timeslice=timeslice)
			else:
				mset.midplaner(selector,skip=skip,rounder=rounder,framecount=framecount,
					protein_selection=protein_select,timeslice=timeslice)
			mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
			mset.analyze_undulations()
			#---Save the data
			pickledump(mset,'pkl.structures.'+tests[testno]+'.'+basename[:11]+'.'+str(timeslice[0])+'-'+
				str(timeslice[1])+'-'+str(timeslice[2])+'.pkl',directory=pickles)
			if erase_when_finished:
				del mset
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


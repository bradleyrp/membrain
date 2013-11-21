#!/usr/bin/python -i

# when running on dirac:
# (1) disable "-i" above and use "-log" 
# (2) comment-out the plotter.py below
# (3) change location flags

from membrainrunner import *
execfile('plotter.py')

#---Set locations, parameters
if 0: # on dirac
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
if 1: # on light
	basedir = '/home/rpb/worker-big/membrane-test-set'
	locations = '/home/rpb/worker/membrain/trajectory-rpb-light'	
if 0: # on dark
	basedir = '/'
	locations = '/store-delta/worker/membrain/trajectory-rpb-dark'

#---Set selections
if 0:
	tests = ['membrane-v509','membrane-v510','membrane-v511']
	ionnames = ['NA','MG','Cal']
	director = ['name P and not resname CHL1','name C218','name C318']
	selector = 'name P'
	residues = ['DOPC','DOPS','PIPP']	
if 1:
	tests = ['membrane-v509','membrane-v510']
	ionnames = ['NA','MG']
	director = ['name P and not resname CHL1','name C218','name C318']
	selector = 'name P'
	residues = ['DOPC','DOPS','PIPP']
if 0:
	tests = ['membrane-v509']
	ionnames = ['NA']
	residues = ['DOPC','DOPS','PI2P']
	director = ['name P and not resname CHL1','name C218','name C318']
	selector = 'name P'
if 0:
	tests = ['membrane-v530','membrane-v531','membrane-v532']
	ionnames = ['NA','MG','Cal']
	residues = ['POPC','CHL1','DOPE','DOPS','PIP2']
	director = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
	selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

skip = 1

#---Analyze the data by looping over systems and trajectories.

starttime = time.time()
[systems,structures,trajectories] = parse_locations_file(basedir,locations)	
print 'Starting analysis job.'
for t in range(len(tests)):
	print 'Running calculation: average structure and undulation spectrum on '+tests[t]+'.'
	for traj in trajectories[systems.index(tests[t])][-1:]:
		#---Load the trajectory
		gro = structures[systems.index(tests[t])]
		basename = traj.split('/')[-1][:-4]
		sel_surfacer = sel_aamd_surfacer
		print 'Accessing '+basename+'.'
		mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
		#---Average structure calculation
		mset.identify_monolayers(director,startframeno=3)
		mset.midplaner(selector,skip=skip,rounder=4.0)
		mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
		mset.analyze_undulations(redundant=0)
		#---Save the data
		pickledump(mset,'pkl.avgstruct.'+tests[t]+'.'+basename+'.pkl')
		#---Clear the class object
		#del mset
		#mset = MembraneSet()
print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'

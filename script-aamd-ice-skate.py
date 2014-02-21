#!/usr/bin/python

# Pseudocode.
# Needs: surface pickle + ion trajectory
# 1. For each ion, from start to end by time (triple loop)
# 		Record the distance from start point by interval length
#		Record if in the zone in separate array (1 or 0)
# 2. Filter array (1 or 0) by some % time in the zone
# 		Only include delta t distances if the mean in the 1/0 array is above some level

interact = True
from membrainrunner import *
execfile('locations.py')

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = None

#---parameters
rounder = 4.0

#---selections
director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector_aamd_symmetric = 'name P'
selector_aamd_asymmetric = '(name P and not resname CHL1)'
selector_aamd_asymmetric = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
residues_aamd_symmetric = ['DOPC','DOPS','PI2P']
residues_aamd_asymmetric = ['DOPC','DOPS','DOPE','POPC','P35P','PI2P']
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']
cgmd_protein = 'name BB'

#---possible analyses
analysis_descriptors = {
	'v510-40000-90000-100':
		{'sysname':'membrane-v510',
		'sysname_lookup':'membrane-v510-atomP',
		'director':director_aamd_symmetric,'selector':selector_aamd_symmetric,'protein_select':None,
		'trajsel':'s8-kraken-md.part0021.40000-90000-100.atomP.xtc',
		'timeslice':[40000,90000,100],
		'ions_struct',
		'ions_traj'}}
analysis_names = ['v510-40000-90000-100']

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def trajectory_lookup(aname):
	'''Return the correct trajectory and structure files.'''
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	if 'sysname_lookup' in vars() and sysname_lookup == None: sysname_lookup = sysname
	grofile = structures[systems.index(sysname_lookup)]
	if type(trajsel) == slice:
		trajfile = trajectories[systems.index(sysname_lookup)][trajsel]
	elif type(trajsel) == str:
		pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajsel))
		for fname in trajectories[systems.index(sysname_lookup)]:
			if pat.match(fname):
				trajfile = [fname]
	elif type(trajsel) == list:
		trajfile = []
		for trajfilename in trajsel:
			pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajfilename))
			for fname in trajectories[systems.index(sysname_lookup)]:
				if pat.match(fname):
					trajfile.append(fname)
	return grofile,trajfile

#---MAIN
#-------------------------------------------------------------------------------------------------------------

starttime = time.time(); print 'start'
#---loop over analysis questions
for aname in analysis_names:
	grofile,trajfile
	#---loop over trajectory files
	for traj in trajfile:
		mset = MembraneSet()
		#---Load the trajectory

		basename = traj.split('/')[-1][:-4]
		#---revised basename to include step-part because sometimes the time gets reset
		basename = "-".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])[:-4]
		sel_surfacer = sel_aamd_surfacer
		print 'accessing '+basename+'.'
		starttime = time.time()
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		print 'time = '+str(1./60*(time.time()-starttime))+' minutes.'

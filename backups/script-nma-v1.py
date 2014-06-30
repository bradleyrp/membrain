#!/usr/bin/python -i

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
location = 'light'
skip = None
framecount = 30

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-repository'
	locations = '/home/rpb/worker/membrain/trajectory-rpb-light'	
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/'
	locations = '/store-delta/worker/membrain/trajectory-rpb-dark'
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)	

if 0:
	from membrainrunner import *
	sysno = systems.index('membrane-v614')
	trajno = -1
	traj = trajectories[sysno][trajno]
	gro = structures[sysno]
	basename = traj.split('/')[-1][:-4]
	mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='cgmd')

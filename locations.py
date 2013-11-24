#!/usr/bin/python 

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
elif location == 'light':
	basedir = '/home/rpb/worker-big/membrane-repository/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-light'	
	pickles = '/home/rpb/worker-big/membrane-repository/pickle-repository/'
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/'
	locations = '/store-delta/worker/membrain/locations-rpb-trajectory-dark'
	pickles = '/store-delta/worker/repo-pickle/'
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

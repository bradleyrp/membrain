#!/usr/bin/python 

#---Location-specific settings
if location == 'dirac':
	basedir = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/'
	locations = '/home/ryb/storedrive/ryb/membrane-v5xx-master-copy/trajectory-map-membrane-v5xx'
	erase_when_finished = True
	plot_suppress = True
elif location == 'light':
	basedir = '/home/rpb/worker/repo-membrane/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-light'	
	pickles = '/home/rpb/worker/repo-pickles/'
	plot_suppress = False
	execfile('plotter.py')
	erase_when_finished = False
elif location == 'dark':
	basedir = '/'
	locations = '/store-delta/worker/membrain/locations-rpb-trajectory-dark'
	pickles = '/store-delta/worker/repo-pickle/'
	plot_suppress = False
	execfile('plotter.py')
	erase_when_finished = False
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

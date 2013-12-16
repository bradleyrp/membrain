#!/usr/bin/python 

#---Automatically detect location
if location == '':
	if subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'dark.site': location = 'dark'
	elif subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'light.site': location = 'light'
	elif subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'dirac': location = 'dirac'

#---Set data locations according to system: dirac
if location == 'dirac':
	basedir = '/media/store-pascal/ryb/membrane-v5xx/'
	locations = '/media/store-pascal/ryb/membrane-v5xx/trajectory-map-membrane-v5xx'
	pickles = '/media/store-pascal/ryb/worker/repo-pickles/'
	erase_when_finished = True
	plot_suppress = True
#---Set data locations according to system: RPB laptop
elif location == 'light':
	basedir = '/home/rpb/worker/repo-membrane/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-light'	
	pickles = '/home/rpb/worker/repo-pickles/'
	plot_suppress = False
	execfile('plotter.py')
	erase_when_finished = False
#---Set data locations according to system: RPB desktop
elif location == 'dark':
	basedir = '/store-delta/compbio/'
	locations = '/store-delta/worker/membrain/locations-rpb-trajectory-dark'
	pickles = '/store-delta/worker/repo-pickles/'
	plot_suppress = False
	execfile('plotter.py')
	erase_when_finished = False

#---Load locations from the table-of-contents
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

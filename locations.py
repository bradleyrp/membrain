#!/usr/bin/python 

#---match hostname with a location script in the locations folder
location_scripts = {
	'dark.site':'./locations/locations-settings-rpb-dark',
	'light.site':'./locations/locations-settings-rpb-light',
	'dirac':'./locations/locations-settings-dirac',
	'arctic':'./locations/locations-settings-ds',
	}

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---execute a location-specific settings file according to hostname
#---the location-specific settings file will load key variables, and also load a directory of simulations
not_barebones = False
if 'location' not in globals() or location == '' or location == None:
	not_barebones = False
	for key in location_scripts.keys():
		# if subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == key:
		# This does not always work when running in certain terminals. Stack says this is better. http://goo.gl/tcNpKc
		import socket
		if socket.gethostname() == key:
			print('status: Loading site-specific ' + key + ' location files.')
			execfile(location_scripts[key])
			not_barebones = True
			break

#---if no locations-settings file is specified, use default parameters
if not not_barebones:
	status('status: using default locations in the current working directory (see locations.py)')
	basedir = os.path.abspath(os.path.expanduser('.'))
	locations = os.path.abspath(os.path.expanduser('./'))
	pickles = os.path.abspath(os.path.expanduser('.'))

#---load locations from the table-of-contents from the parse function defined in membrainrunner.py

# It seems basedir is only defined it not_barebones is False (lines 28, 30)
if not_barebones: [systems,structures,trajectories] = parse_locations_file(basedir,locations)

#---register the postmortem function in membrainrunner.py with the list of globals
#---whenever interact is set or sent as a flag (-i) the program will conclude with an interactive terminal
if 'interact' in globals() and interact and interact_registered == False:
	atexit.register(postmortem,banner='status: here is an interactive terminal',scriptglobals=globals())
	interact_registered = True

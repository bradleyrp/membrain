#!/usr/bin/python -i

#---match hostname with a location script in the locations folder
#---if this is executed from a directory not called "membrain", then it must be a subdirectory
if '__file__' in globals() and os.path.dirname(os.path.realpath(__file__)).split('/')[-1] != 'membrain':
	dir_prefix = '../'
else: dir_prefix = './'
location_scripts = {
	'dark.site':dir_prefix+'locations/locations-settings-rpb-dark',
	'light.site':dir_prefix+'locations/locations-settings-rpb-light',
	'dirac':dir_prefix+'locations/locations-settings-dirac',
	'ds':dir_prefix+'locations/locations-settings-ds',
	}

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---execute a location-specific settings file according to hostname
#---the location-specific settings file will load key variables, and also load a directory of simulations
not_barebones = False
if 'location' not in globals() or location == '' or location == None:
	not_barebones = False
	for key in location_scripts.keys():
		if subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == key:
			execfile(location_scripts[key])
			not_barebones = True

#---if no locations-settings file is specified, use default parameters
if not not_barebones:
	status('status: using default locations in the current working directory (see locations.py)')
	basedir = os.path.abspath(os.path.expanduser('.'))
	locations = os.path.abspath(os.path.expanduser('./'))
	pickles = os.path.abspath(os.path.expanduser('.'))

#---load locations from the table-of-contents from the parse function defined in membrainrunner.py
if not_barebones: [systems,structures,trajectories] = parse_locations_file(basedir,locations)

#---register the postmortem function in membrainrunner.py with the list of globals
#---whenever interact is set or sent as a flag (-i) the program will conclude with an interactive terminal
if 'interact' in globals() and interact and 'interact_registered' not in globals():
	#---note that I previously tried to make interact_registered explicit to prevent multiple registrations
	#---...this seemed to fail because importing membrainrunner.py did not notice that interact_registered
	#---...had been previously loaded into globals and then set to true, so there was clearly a scoping issue
	#---...to get around this I just use the presence of interact_registered in globals as a sign that one
	#---...registration has occured however this is admittedly a bit much for a simple terminal function
	atexit.register(postmortem,banner='status: here is an interactive terminal',scriptglobals=globals())
	interact_registered = True

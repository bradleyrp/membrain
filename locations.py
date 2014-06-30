#!/usr/bin/python 

#---automatically detect location
if 'location' not in globals() or location == '' or location == None:
	if subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'dark.site':
		location = 'dark'
	elif subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'light.site':
		location = 'light'
	elif subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'dirac':
		location = 'dirac'
	elif subprocess.check_output(['echo $HOSTNAME'],shell=True).strip('\n') == 'ground-control':
		location = 'ds'

#---dirac, location specific settings
if location == 'dirac':
	basedir = '/media/store-pascal/ryb/membrane-v5xx/'
	locations = '/media/store-pascal/ryb/membrane-v5xx/trajectory-map-membrane-v5xx'
	pickles = '/media/store-pascal/ryb/worker/repo-pickles/'
	erase_when_finished = True
	plot_suppress = True
#---rpb (laptop), location specific settings
elif location == 'light':
	basedir = '/home/rpb/worker/repo-membrane/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-light'	
	pickles = '/home/rpb/worker/repo-pickles/'
	plot_suppress = False
	if os.path.exists('plotter.py'): execfile('plotter.py')
	else: import matplotlib as mpl
	#---commands for sans-serif fonts on all plots
	mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	mpl.rc('text', usetex=True)
	if 0:
		mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	else:
		#---hacks to allow "boldsymbol" on wispy kappas
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{sfmath}',r'\usepackage{amsmath}',r'\usepackage{siunitx}',r'\sisetup{detect-all}',
			r'\usepackage{helvet}',r'\usepackage{sansmath}',r'\sansmath']  
	erase_when_finished = False
	fsaxtext = 18
	fsaxlabel = 18
	fsaxticks = 18
	fsaxtitle = 18
#---rpb (desktop), location specific settings
elif location == 'dark':
	basedir = '/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-dark'
	pickles = '/home/rpb/worker/repo-pickles/'
	plot_suppress = False
	if os.path.exists('plotter.py'): execfile('plotter.py')
	else: import matplotlib as mpl
	#---commands for sans-serif fonts on all plots
	mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	mpl.rc('text', usetex=True)
	if 0: mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	else:
		#---hacks to allow "boldsymbol" on wispy kappas
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{sfmath}',r'\usepackage{amsmath}',r'\usepackage{siunitx}',r'\sisetup{detect-all}',
			r'\usepackage{helvet}',r'\usepackage{sansmath}',r'\sansmath']  
	erase_when_finished = False
	fsaxtext = 16
	fsaxlabel = 16
	fsaxticks = 16
	fsaxtitle = 20
	fsaxlegend = 14
	fsaxlegend_small = 12
	#---distinct from plot_suppress rpb added this for remote plotting on a specific set of scripts
	plotviewflag = True
	showplots = False
#---ds, location specific settings
elif location == 'ds':
	basedir = '/home/davids/membrane-v5xx/'
	locations = '/home/davids/membrain/locations-ds'
	pickles = '/home/davids/repo-pickles/'
	plot_suppress = False
	execfile('plotter.py')
	#---plot commands previously found in plotter.py
	mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	mpl.rc('text', usetex=True) 
	mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	fsaxtext = 16
	fsaxlabel = 16
	fsaxticks = 16
	fsaxtitle = 20
	fsaxlegend = 14
	if 0: font = {'family' : 'sans-serif', 'size'   : 22}
	font = {'family' : 'sans-serif'}
	mpl.rc('font', **font)
	mpl.rc('text', usetex=True)
	mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	mpl.rcParams['text.latex.preamble'] = [
		r'\usepackage{sfmath}',
		r'\usepackage{amsmath}',
		r'\usepackage{siunitx}',
		r'\sisetup{detect-all}',
		r'\usepackage{helvet}',
		r'\usepackage{sansmath}',
		r'\sansmath',
		r'\usepackage{upgreek}']
	if 0: mpl.rcParams['xtick.major.pad'] = 8
	if 0: mpl.rcParams['ytick.major.pad'] = 8

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load locations from the table-of-contents from the parse function defined in membrainrunner.py
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

#---register the postmortem function in membrainrunner.py with the list of globals
#---whenever interact is set or sent as a flag (-i) the program will conclude with an interactive terminal
if 'interact' in globals() and interact:
	atexit.register(postmortem,
		banner='status: interactive terminal',
		scriptglobals=globals())


#!/usr/bin/python 

#---Automatically detect location
if 'location' not in globals() or location == '' or location == None:
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
#---Set data locations according to system: RPB desktop
elif location == 'dark':
	basedir = '/'
	locations = '/home/rpb/worker/membrain/locations-rpb-trajectory-dark'
	pickles = '/home/rpb/worker/repo-pickles/'
	plot_suppress = False
	execfile('plotter.py')
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
#---Set data locations according to system: DS
elif location == 'ds':
	#---Nb: put system-specific commands here. 
	basedir = ''
	locations = ''
	pickles = ''
	plot_suppress = False
	execfile('plotter.py')
	#---plot commands previously found in plotter.py
	mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	mpl.rc('text', usetex=True) 
	mpl.rc('text.latex', preamble='\usepackage{sfmath}')

#---Load locations from the table-of-contents
[systems,structures,trajectories] = parse_locations_file(basedir,locations)

#---Post-mortem (cleanup) function
#---Nb: as far as I can tell, there is no way to implicitly send the script globals() object back to
#---our wrapper module (membrainrunner.py). The other direction is easy with either "from __main__ import *" or 
#---simply using execfile, in which case an import/execfile command drops the globals into the module. The other
#---direction is impossible, and since we always call locations.py with execfile, it makes sense to put the
#---following interactive terminal option in here. When you register the postmortem function here, it doesn't
#---actually take globals until the script tries to exit, so it doesn't matter that this comes early.
#---Would be nice to find a way to use the exception trick to grab the namespace, or something similar.
if 'interact' in globals() and interact:
	if 'postmortem' not in [i[0].__name__ for i in atexit._exithandlers]:
		atexit.register(postmortem,
			banner='Here is the interactive terminal you wanted.',
			scriptglobals=globals())

#---Universal definitions
#---Nb: these are used for parsing frames, but I set the defaults here
skip = None
framecount = None


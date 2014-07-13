#!/usr/bin/python

#---MEMBRAINRUNNER: A Convenient wrapper for the MEMBRAIN bilayer analysis package.

"""
Needs better documentation.
Possibly needs more generic name like "simcubate".
"""

#---LIBRARIES
#-------------------------------------------------------------------------------------------------------------

#---note: the following import is unclear. it may have been used for a different module import arrangement
from __main__ import *
import sys, atexit, code
#---debug mode can be enabled from the command line
if '-d' in sys.argv: debugmode = True
#---interactive mode can be enabled from the command line
if '-i' in sys.argv: interact = True
#---override the logfile if interactive mode is enabled
if ('-i' in sys.argv or ('interact' in globals() and interact)) \
	and 'logfile' in globals() and logfile != None: 
	logfile = None
interact_registered = False

#---import the primary membrain library for which membrainrunner is only a wrapper
from membrain import *

#---declare an empty MembraneSet object if none exists
#---temporarily disabled since this is redundant
if 0 and 'mset' not in globals(): mset = MembraneSet()

#---INTERACTIVE MODE HANDLERS
#-------------------------------------------------------------------------------------------------------------

class ErrorHandler:
	'''Deprecated function which logged errors if necessary, but was not working with streams properly.'''
	def __init__(self,errlogname):
		self.errlogname = errlogname
	def write(self,string):
		fname = self.errlogname
		handler = open(fname,'w',0)
		handler.write(string)
		handler.close()
	def close(self): pass
	def flush(self): pass
		
def test(args):
	'''Test function.'''
	print '\n\nThis is a test! Watch out for occult enemies !!!\n'
	print 'Previously, we used wrappers in membrain.sh and membrain-tool.sh to call functions in ',
	print 'membrainrunner.py. These days, we just write self-contained scripts, but you can still find ',
	print 'the tools in the beta folder. They make for a nice constellation of scripts for integrating ',
	print 'python wrappers, bash wrappers, nice compact bash flags, an optional logging system, ',
	print 'and an interface which pulls up data files depending on which system you are using. ',
	print 'And so on. \n\n\t-Ryan Bradley 2014.02.13\n'
	
def postmortem_debug(type,value,tb):
	#---if already in interactive mode or no tty-like device call default hook
	if hasattr(sys, 'ps1') or not sys.stderr.isatty():
		sys.__excepthook__(type, value, tb)
	#---not in interactive mode
	else:
		import traceback,pdb
		traceback.print_exception(type,value,tb)
		print
		pdb.pm()

def postmortem(scriptglobals=None,banner=None):
	'''Post-mortem clean-up function which supplies a prompt via atexit and locations.py.
	
	This supplies an interactive prompt.
	It receives the globals via the atexit.register function.
	Seemingly impossible to pass globals back to an imported module.
	So instead we just call atexit.register(postmortem,scriptglobals=globals()) from locations.py.
	Which is always called via execfile in the script. 
	This trick works well enough.
	'''
	code.interact(local=scriptglobals,banner=banner)

#---TRAJECTORY LIBRARY HANDLERS
#-------------------------------------------------------------------------------------------------------------

"""
The following funtions are designed to easily lookup molecular dynamics trajectories from a table-of-contents.
The parse_locations_file is called in locations.py which is called via execfile to load all the paths.
"""
		
def parse_locations_file(basedir,locations):
	'''Parse the trajectory location file.'''
	systems = []
	trajectories = []
	structures = []
	fp  = open(locations,'r')
	line = fp.readline()
	while line:
		location_set = []
		if line[0:5] == 'name=':
			systems.append(line.strip('\n')[5:])
			line = fp.readline()
			location_set = []
			while line:
				if line[0:7] == 'struct=':
					structures.append(line[7:].strip())
				elif line[0:5] == 'name=':
					trajectories.append(location_set)
					location_set = []
					break
				elif line[0] != '#' and line[0].strip('\n') != '':
					location_set.append(line.strip())
				line = fp.readline()
		else: line = fp.readline()
		if location_set != []: 
			trajectories.append(location_set)
			location_set = []
	fp.close()
	return [systems,structures,trajectories]
	
def trajectory_lookup(analysis_descriptors,aname,globs,
	keytrajsel='trajsel',keysysname='sysname_lookup'):
	'''Return the correct trajectory and structure files, given a dictionary and key.'''
	#---Since we are in an imported module (membrainrunner.py), we have to get globals for the lookup
	structures = globs['structures']
	systems = globs['systems']
	trajectories = globs['trajectories']
	#---Also lookup entries from the dictionary
	if keysysname not in analysis_descriptors[aname].keys(): sysname_lookup = None
	else: sysname_lookup = (analysis_descriptors[aname])[keysysname]
	trajsel = (analysis_descriptors[aname])[keytrajsel]
	#---If sysname isn't defined, lookup the proper sysname from the dictionary
	if sysname_lookup == None: sysname_lookup = (analysis_descriptors[aname])['sysname']
	grofile = structures[systems.index(sysname_lookup)]
	#---Nb this routine searches the locations file for the right trajectory file.
	#---Requires analysis descriptors dictionary of dictionaries with the following key variables
	#--->	1. sysname_lookup : key for the system name, i.e. membrane-v510-atomP
	#--->	2. trajsel : a variable type which specifies the trajectory file you want (see below)(
	if type(trajsel) == slice:
		trajfile = trajectories[systems.index(sysname_lookup)][trajsel]
	elif type(trajsel) == str:
		pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajsel))
		for fname in trajectories[systems.index(sysname_lookup)]:
			if pat.match(fname): trajfile = [fname]
	elif type(trajsel) == list:
		trajfile = []
		for trajfilename in trajsel:
			pat = re.compile('(.+)'+re.sub(r"/",r"[/-]",trajfilename))
			for fname in trajectories[systems.index(sysname_lookup)]:
				if pat.match(fname):
					trajfile.append(fname)
	else:
		trajfile = None
		grofile = None
	return grofile,trajfile
	
def specname_from_pickle(basename):
	'''Standard method for deriving a step/part/timeslice-specific name from the trajectory file.'''
	special_name = ".".join([i for i in re.match('.+membrane\-.+',basename).string.split('.')
		if re.match('[a-z][0-9]\-.+',i) 
		or re.match('membrane-v.+',i) 
		or re.match('md',i)
		or re.match('part[0-9]{4}',i)
		or re.match('[0-9]+\-[0-9]+\-[0-9]',i)])
	return special_name

def specname_from_traj(traj):
	'''Return basic name from a trajectory file name.'''
	#---Nb previous version might include subset names (like 'atomP') but now we drop them
	#return ".".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])[:-4]
	return '.'.join((".".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])).split('.')\
		[:(-2 if len(traj.split('.')[-2].split('-')) == 1 else -1)])

def specname_guess(sysname,trajsel):
	'''If you have a trajsel but it's not in the table of contents you can guess the specname.'''
	#---Nb previous version might include subset names (like 'atomP') but now we drop them
	#return ".".join(re.match('.*/[a-z][0-9]\-.+',traj).string.split('/')[-2:])[:-4]
	return sysname+'.'+'.'.join((".".join(re.match('[a-z][0-9]\-.+',trajsel).string.split('/')\
		[-2:])).split('.')[:(-2 if len(trajsel.split('.')[-2].split('-')) == 1 else -1)])

def specname_pickle(sysname,traj,timeslice=None):
	'''Construct a standard picklename from a system name name and trajectory name.'''
	if timeslice == None: picklename = sysname+'.'+specname_from_traj(traj)
	else: picklename = sysname+'.'+('.'.join(specname_from_traj(traj).split('.')[:-1]))+\
		('-'.join([str(i) for i in timeslice]))
	return picklename

#---MONITORING FUNCTIONS
#-------------------------------------------------------------------------------------------------------------
	
def status(string,start=None,i=None,looplen=None):
	'''Print status to the screen also allows for re-writing the line. Duplicated in the membrain library.'''
	#---note: still looking for a way to use the carriage return for dynamic counter without ...
	#---...having many newlines printed to the file. it seems impossible to use the '\r' + flush()...
	#---...method with both screen and file output, since I can't stop the buffer from being written
	#---display a refreshable string
	if start == None and looplen == None and i != None:		
		print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8)+'\n',
	elif start == None and looplen != None and i != None:		
		if i+1 == looplen:
			print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8)+'\n',
		#---if a logfile has been defined, this output is destined for a file in which case suppress counts
		elif i+1 != looplen and ('logfile' not in globals() or logfile == None):
			print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8),
			sys.stdout.flush()
	#---estimate the remaining time given a start time, loop length, and iterator
	elif start != None and i != None and looplen != None and ('logfile' not in globals() or logfile == None):
		esttime = (time.time()-start)/(float(i+1)/looplen)
		print '\r'+string.ljust(20)+str(abs(round((esttime-(time.time()-start))/60.,1))).ljust(10)+\
			'minutes remain',
		sys.stdout.flush()
	#---standard output
	else: print string
	
#---record a global start time
global_start_time = time.time()

def checktime():
	'''Report the current time.'''
	status('status: time = '+str(1./60*(time.time()-global_start_time))+' minutes')

#---MAIN
#-------------------------------------------------------------------------------------------------------------

'''
This section handles incoming arguments, logging, error logging, and inducing a debug mode if desired.
'''

#---completely override argument handling if it already exists
if 'args' not in globals():
	#---parser and logging
	parser = argparse.ArgumentParser(description='Membrain argument parser.',prog='Membrain')
	parser.add_argument('-viz',action='store_true',help='Set this to use mayavi 3D visualization.')
	parser.add_argument('-d','--debugmode',action='store_true',help='Set this to use mayavi 3D visualization.')
	parser.add_argument('-log',help='The log file. Standard output routes here. Error goes to a separate file.')
	parser.add_argument('-o','--operation',choices=['test'],help='Calculation to perform.')
	parser.add_argument('-i','--interactive',action='store_true',help='Run python in interactive mode.')
	args = parser.parse_args()

	#---enable 3D visualization
	if args.viz == 'yes' or args.viz == True:
		from mayavi import mlab

	#---write stdout to a log 
	#---enable via the "-log" flag to either the script or membrainrunner
	#---enable also via setting the logfile variable to a string in the script
	if args.log != None or ('logfile' in globals() and logfile != None):
		import sys
		stdstderr = sys.stderr
		if args.log != None: logfile = args.log
		sys.stdout = sys.stderr = open(logfile,'w',1)

	#---switching to interactive mode on error
	#---enable via the "-d" flag to either the script or membrainrunner
	#---enable also via setting the debugmode variable to True in the script
	if ('debugmode' in globals() and debugmode) or args.debugmode:
		sys.excepthook = postmortem_debug
	
	#---call any function called when the operation flag is caught by the parser
	if __name__ == "__main__" and args.operation != None:
		globals()[args.operation](args)


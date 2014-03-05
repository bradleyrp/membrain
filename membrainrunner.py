#!/usr/bin/python

#---MEMBRAINRUNNER: A Convenient wrapper for the MEMBRAIN bilayer analysis package.

#---LIBRARIES
#-------------------------------------------------------------------------------------------------------------

from __main__ import *
import sys, atexit, code
if '-d' in sys.argv: debugmode = True
if '-i' in sys.argv:
	interact = True
	if 'logfile' in globals() and logfile != None:
		print 'overriding log file because you set interactive flag'
		logfile = None

from membrain import *

if 'mset' not in globals(): mset = MembraneSet()

nspass = None
#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

class ErrorHandler:
	def __init__(self,errlogname):
		self.errlogname = errlogname
	def write(self,string):
		fname = self.errlogname
		handler = open(fname,'w')
		handler.write(string)
		handler.close()
	def close(self):
		pass
	def flush(self):
		pass
		
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
	
def trajectory_lookup(analysis_descriptors,aname,globs):
	'''Return the correct trajectory and structure files, given a dictionary and key.'''
	#---Since we are in an imported module (membrainrunner.py), we have to get globals for the lookup
	structures = globs['structures']
	systems = globs['systems']
	trajectories = globs['trajectories']
	#---Also lookup entries from the dictionary
	if 'sysname_lookup' not in analysis_descriptors[aname].keys(): sysname_lookup = None
	else: sysname_lookup = (analysis_descriptors[aname])['sysname_lookup']
	trajsel = (analysis_descriptors[aname])['trajsel']
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
    '''Post-mortem clean-up function which supplies a prompt via atexit and locations.py.'''
	#---This supplies an interactive prompt.
	#---It receives the globals via the atexit.register function.
	#---Seemingly impossible to pass globals back to an imported module.
	#---So instead we just call atexit.register(postmortem,scriptglobals=globals()) from locations.py.
	#---Which is always called via execfile in the script. This trick works well enough.
    code.interact(local=scriptglobals,banner=banner)

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
	'''Construct a standard picklename from a systename name and trajectory name.'''
	if timeslice == None:
		picklename = sysname+'.'+specname_from_traj(traj)
	else:
		picklename = sysname+'.'+('.'.join(specname_from_traj(traj).split('.')[:-1]))+\
			('-'.join([str(i) for i in timeslice]))
	return picklename

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---parser and logging
parser = argparse.ArgumentParser(description='Membrain argument parser.',prog='Membrain')
parser.add_argument('-viz',action='store_true',help='Set this to use mayavi 3D visualization.')
parser.add_argument('-d','--debugmode',action='store_true',help='Set this to use mayavi 3D visualization.')
parser.add_argument('-log',help='The log file. Standard output routes here. Error goes to a separate file.')
parser.add_argument('-o','--operation',choices=['test'],help='Calculation to perform.')
parser.add_argument('-i','--interactive',action='store_true',help='Run python in interactive mode.')
args = parser.parse_args()

#---enable visualization
if args.viz == 'yes' or args.viz == True:
	from plotter import *

#---write stdout to a log 
#---enable via the "-log" flag to either the script or membrainrunner
#---enable also via setting the logfile variable to a string in the script
if args.log != None:
	import sys
	stdstderr = sys.stderr
	sys.stderr = ErrorHandler(args.log+'.err')
	sys.stdout = open(args.log,'w',0)
elif 'logfile' in globals():
	if logfile != None:
		import sys
		stdstderr = sys.stderr
		sys.stderr = ErrorHandler(logfile+'.err')
		sys.stdout = open(logfile,'w',0)

#---switching to interactive mode on error
#---enable via the "-d" flag to either the script or membrainrunner
#---enable also via setting the debugmode variable to True in the script
if ('debugmode' in globals() and debugmode) or args.debugmode:
	sys.excepthook = postmortem_debug
	
#---global start time
global_start_time = time.time()

#---report the current time
def checktime():
	print 'status: time = '+str(1./60*(time.time()-global_start_time))+' minutes'

#---pass arguments and call a function
if __name__ == "__main__" and args.operation != None:
	globals()[args.operation](args)


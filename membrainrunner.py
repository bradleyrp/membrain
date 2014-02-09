#!/usr/bin/env python

#---LIBRARIES
#-------------------------------------------------------------------------------------------------------------

from membrain import *
execfile('/etc/pythonstart')

#---VARIABLES
#-------------------------------------------------------------------------------------------------------------

#---Selection groups, AAMD
lipidres = ['DOPC','DOPS','PI2P']
head1 = 'N  C12  C13  C14  C15 H12A H12B H13A H13B H13C '+\
	'H14A H14B H14C H15A H15B H15C  C11 H11A H11B  P  '+\
	'O13  O14  O12  O11   C1   HA   HB   C2   HS'
head2 = 'N   HN1   HN2   HN3   C12  H12A   C13  O13A  O13B'+\
	'   C11  H11A  H11B     P   O13   O14   O12   '+\
	'O11    C1    HA    HB    C2    HS'
head3 = 'C12   H2   O2  HO2  C13   H3   O3  HO3  C14   H4'+\
	'   O4   P4 OP42 OP43 OP44  C15   H5   O5   P5 '+\
	'OP52 OP53 OP54  H52  C16   H6   O6  HO6  C11   H1'+\
	'    P  O13  O14  O12  O11   C1   HA   HB   C2   HS'
lipidresheads = [head1.split(),head2.split(),head3.split()]
sel_aamd_lipid_heads = []
for j in range(len(lipidres)):
	temp = 'resname '+lipidres[j]+' and (name '
	for i in range(len(lipidresheads[j])-1):
		temp=temp+lipidresheads[j][i]+' or name '
	temp=temp+lipidresheads[j][i+1]+')'
	sel_aamd_lipid_heads.append(temp)
sel_aamd_heads = '('+sel_aamd_lipid_heads[0]+') or ('+sel_aamd_lipid_heads[1]+') or ('+\
	sel_aamd_lipid_heads[2]+')'
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']

#---Selection groups, CGMD
sel_cgmd_surfacer = ['name PO4 or name POG','name C2A']

#---Selection groups, all
sel_all = ['all']

#---Master membrane class object
if 'mset' not in globals(): mset = MembraneSet()

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def undulations(args):
	if args.resolution == 'meso':
		#---Rectangular mesoscale data
		if args.xyzform == 'rect':
			mset.load_points(args.dir,
				xyzform=args.xyzform,
				nbase=int(args.nbase),
				rounder=1.0,
				lenscale=(float(args.length) if args.length != None else None),
				start=(int(args.start) if args.start != None else None),
				skip=(int(args.skip) if args.skip != None else None))
			mset.surfacer(rounder=2.0)
			mset.calculate_undulation_spectrum(removeavg=0,redundant=1)
			mset.analyze_undulations(redundant=1)
			plotter_undulations_summary(mset,qmagfilter=[0.1,1.],zoom=True)
		#---Square mesoscale data
		elif args.xyzform == 'square':
			mset.load_points(args.dir,
				xyzform=args.xyzform,
				nbase=int(args.nbase),
				lenscale=(float(args.length) if args.length != None else None),
				start=(int(args.start) if args.start != None else None),
				shifter=2*float(args.rounder),
				rounder=float(args.rounder),
				prefix=args.prefix,
				suffix=args.suffix)
			mset.surfacer(rounder=float(args.rounder))
			mset.calculate_undulation_spectrum(removeavg=1,redundant=2)
			mset.analyze_undulations(redundant=2)
			plotter_undulations_summary(mset,qmagfilter=[0.03,0.1],zoom=True)
		#---Square mesoscale data, version 2
		elif args.xyzform == 'square2':
			mset.load_points(args.dir,
				xyzform=args.xyzform,
				nbase=int(args.nbase),
				lenscale=(float(args.length) if args.length != None else None),
				start=(int(args.start) if args.start != None else None),
				shifter=2*float(args.rounder),
				rounder=float(args.rounder),
				prefix=args.prefix,
				suffix=args.suffix)
			mset.surfacer(rounder=float(args.rounder))
			mset.calculate_undulation_spectrum(removeavg=1,redundant=2)
			mset.analyze_undulations(redundant=2)
			plotter_undulations_summary(mset,zoom=True)
		#---Precomputed data
		elif args.xyzform == 'regular':
			mset.load_points(args.dir,
				regular=True,
				nbase=int(args.nbase),
				rounder=float(args.rounder),
				lenscale=(float(args.length) if args.length != None else None),
				start=(int(args.start) if args.start != None else None),
				prefix=args.prefix,
				suffix=args.suffix)
			mset.calculate_undulation_spectrum(removeavg=0,redundant=1)
			mset.analyze_undulations(redundant=1)
			plotter_undulations_summary(mset,qmagfilter=[0.02,0.1],zoom=True)
	#---All-atom MD trajectory
	elif args.resolution == 'aamd':
		sel_surfacer = sel_aamd_surfacer
		mset.load_trajectory((args.config,args.traj),
			resolution=args.resolution,
			start=(int(args.start) if args.start != None else None),
			end=(int(args.end) if args.end != None else None),
			skip=(int(args.skip) if args.skip != None else None))
		mset.identify_monolayers(['name P','name C218','name C318'])
		mset.midplaner('name P',
			interp=args.mode,
			start=(int(args.start) if args.start != None else None),
			skip=(int(args.skip) if args.skip != None else None),
			end=(int(args.end) if args.end != None else None),
			rounder=4.0)
		mset.midplaner('name P',skip=skipper,rounder=4.0)
		mset.calculate_undulation_spectrum(removeavg=0,redundant=0)
		mset.analyze_undulations(redundant=0)
		plotter_undulations_summary(mset,qmagfilter=[.1,2.],zoom=True)
	#---Coarse-grained MD trajectory
	elif args.resolution == 'cgmd':
		sel_surfacer = sel_cgmd_surfacer
		mset.load_trajectory((args.config,args.traj),
			resolution=args.resolution,
			start=(int(args.start) if args.start != None else None),
			end=(int(args.end) if args.end != None else None),
			skip=(int(args.skip) if args.skip != None else None))
		mset.identify_monolayers(['name PO4','name C2A'])
		mset.midplaner('name PO4',
			interp=args.mode,
			start=(int(args.start) if args.start != None else None),
			skip=(int(args.skip) if args.skip != None else None),
			end=(int(args.end) if args.end != None else None),
			rounder=20.0)
		mset.calculate_undulation_spectrum(redundant=1,removeavg=1)
		mset.analyze_undulations(redundant=1)
		plotter_undulations_summary(mset,qmagfilter=[.05,.7],zoom=True)
	return mset

#---Bookkeeping
	
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
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Membrain argument parser.',prog='Membrain')
#---General flags
parser.add_argument('-o','--operation',choices=['gr','undulations'],help='Calculation to perform.')
parser.add_argument('-mode',help='General-purpose parameter, means something different every time.')
parser.add_argument('-start',help='Start frame.')
parser.add_argument('-end',help='End frame.')
parser.add_argument('-skip',help='Skip frames.')
parser.add_argument('-r','--resolution',choices=['aamd','cgmd','meso'],help='Simulation scale.')
#---Flags for MD trajectories
parser.add_argument('-c','--config',help='Starting structure.')
parser.add_argument('-t','--traj',help='Trajectory file.')
#---Flags for the mesoscale code
parser.add_argument('-nbase',help='Size parameter for incoming mesoscale simulation data.')
parser.add_argument('-length',help='Length scale for the input data.')
parser.add_argument('-rounder',help='Generic rounder function, used differently each time.')
parser.add_argument('-xyzform',choices=['rect','square','square2','regular'],help='Which type of xyz data to input.')
parser.add_argument('-d','-dir','--dir',help='Directory location of xyz files.')
parser.add_argument('-prefix',help='When loading xyz trajectories, this is the file prefix.')
parser.add_argument('-suffix',help='When loading xyz trajectories, this is the file suffix.')
#---Code flags
parser.add_argument('-viz',action='store_true',help='Set this to use mayavi 3D visualization.')
parser.add_argument('-log',help='The log file. Standard output routes here. Error goes to a separate file.')
parser.add_argument('-pickle',help='Provide a filename for the pickle output of the calculation.')
parser.add_argument('-locations',help='Provide a formatted list of trajectory locations.')
#---Parser and logging
args = parser.parse_args()
if args.viz == 'yes' or args.viz == True:
	from plotter import *
if args.log != None:
	import sys
	stdstderr = sys.stderr
	sys.stderr = ErrorHandler(args.log+'.err')
	sys.stdout = open(args.log,'w',0)

if __name__ == "__main__":
	globals()[args.operation](args)

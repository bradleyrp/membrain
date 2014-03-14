#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

analysis_descriptors = {
	'v514-10000-29000-100':
		{'sysname':'membrane-v514',
		'sysname_lookup':'membrane-v514-ions',
		'trajsel':'s3-sim-compbio-md.part0004.10000-29000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v514.a2-surfacer.s3-sim-compbio-md.part0004.10000-29000-100.pkl',
		'ionname':'NA'},
	'v532-20000-58000-100':
		{'sysname':'membrane-v532',
		'sysname_lookup':'membrane-v532-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-58000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',
		'ionname':'Cal'},
	'v531-20000-62000-100':
		{'sysname':'membrane-v531',
		'sysname_lookup':'membrane-v531-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-62000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v531.a6-surfacer.s4-sim-trestles-md.part0007.20000-62000-100.pkl',
		'ionname':'MG'},
	'v530-30000-100000-100':
		{'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-ions',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl',
		'ionname':'NA'},
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl',
		'ionname':'Cal'},
	'v509-30000-90000-100':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509',
		'trajsel':'s4-kraken-md.part0007.30000-90000-100.ions.xtc',
		'ionname':'Na',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl'}}
analysis_names = ['v532-20000-58000-100','v530-30000-100000-100','v531-20000-62000-100',
	'v511-30000-80000-100','v509-30000-90000-100'][-1:]
routine = ['compute','postproc','computexyz',][:1]
plot_suppress = True
#---method parameters
upto = 500 #---how far to only look at the diffusion curves

#---method parameters
upto = 500 #---how far to only look at the diffusion curves
timelimit = 2*10**2 #---maximum ps for fitting which depends on occupancy because displacements shrink
occupancy = 1. #---occupancy rate for consideration of the displacement
bwid = 5 #---angstrom width of slices where default is 10 and zonesub must be none
flush_bin_edges = True #---new flush bin edges method is recommended
if not flush_bin_edges: #---never define zonesub manually of using flush method
	zonesub = [1,2,3,4,5,6,7,8,14,15,16,17,18,19,20,21] #---custom selection of zones
	zonesub = None #---only set this manually if flush_bin_edges is turned off
dtlimit = 350 #---time maximum for computation, above which we don't include the data to save memory
ion_drop_thresh = 0.001 #---threshold ion occupancy for analysis
sideslices = 6 #---number of slices to include on each side

#---plot details
nbins = 40
allalpha = 1.
edgeprop = 2.

bwid = 2

#---MAIN
#-------------------------------------------------------------------------------------------------------------
if plot_suppress:
	mpl.use('Agg')
	
#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice
if 'compute' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---check for pre-existing pickle
		resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		#---disabled file-checking for now
		if 0: 
			mdionskate = unpickle(pickles+resultpkl)
			if mdionskate != None: raise Exception('except: pkl already exists so fix naming')
		#---get ion positions and times
		print 'status: getting ions'
		clock = []
		ionsp = []
		ionsn = []
		vecs = []
		ion_select_p = mset.universe.selectAtoms('name '+'NA')
		ion_select_n = mset.universe.selectAtoms('name '+'CL')
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionsp.append(ion_select_p.coordinates())
			ionsn.append(ion_select_n.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
			vecs.append(mset.universe.dimensions[:3])
		ionsnz = array(ionsn)[:,:,2].T
		ionspz = array(ionsp)[:,:,2].T
		nzones = int(mean(vecs,axis=0)[2]/bwid)
		binws = [linspace(0,i[2],nzones)[1] for i in vecs]
		roundzn = array([[int(ionsnz[i][j]/binws[j]) for j in range(len(ionsnz[i]))] 
			for i in range(len(ionsnz))])
		roundzp = array([[int(ionspz[i][j]/binws[j]) for j in range(len(ionspz[i]))] 
			for i in range(len(ionspz))])
		plt.grid(True)
		plt.hist(array(roundzp).flatten(),range=(0,nzones),bins=nzones,alpha=0.5)
		plt.hist(array(roundzn).flatten(),range=(0,nzones),bins=nzones,alpha=0.5)
		plt.show()


#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice
if 'computex' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#mset_surf = unpickle(pickles+structure_pkl)
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---check for pre-existing pickle
		resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		#---disabled file-checking for now
		if 0: 
			mdionskate = unpickle(pickles+resultpkl)
			if mdionskate != None: raise Exception('except: pkl already exists so fix naming')
		#---get ion positions and times
		print 'status: getting ions'
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		#---check for same timestamps and number of frames
		if mset_surf.surf_time != mset.time_list:
			raise Exception('except: timestamps don\'t match')
		#---Nb ds removed the following due to mismatch problem but timestamps match 2014.03.05
		if 0: ionspos = array(ionspos)[:-1]
		print len(ionspos)
		ionstraj = []
		for ind in range(shape(ionspos)[1]):
			course = array(ionspos)[:,ind]
			#---three-line handling PBCs
			hoplistp = (course[1:]-course[:-1])>array(mset_surf.vecs)[1:]/2.
			hoplistn = (course[:-1]-course[1:])>array(mset_surf.vecs)[1:]/2.
			course_nojump = course[1:]-(cumsum(1*hoplistp-1*hoplistn,axis=0))*array(mset_surf.vecs)[1:]
			ionstraj.append(course_nojump)
		ionstraj=array(ionstraj)
		nions = len(ionstraj)
		nframes = len(ionstraj[0])
		if dtlimit > nframes:
			dtlimit = nframes
			print 'warning: dtlimit was higher than nframes so watch out for statistical anomalies'
		#---specify zones
		center = mean(mset_surf.surf_position)
		cslice = int(round(center/bwid))
		thick = mean(mset_surf.surf_thick)
		zonesabs = list(arange(center,0-bwid,-bwid))[::-1][:-1]+list(arange(center,vecs[2]+bwid,bwid))
		zonecenters = array(zonesabs[:-1]+zonesabs[1:])/2.
		#---handle the bin arithmetic
		if not flush_bin_edges:
			#---original method here, before fixing flush bin edges
			rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+\
				(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)]
			roundz = [[int(round(i)) for i in (rewrapz[j]-center)/bwid] for j in range(nions)]
			roundz = roundz - array(roundz).min()
		else:
			#---new method ensures flush bin edges and sets zonesub
			nzones = int(mean(mset_surf.vecs,axis=0)[2]/bwid)
			binws = [linspace(0,i[2],nzones)[1] for i in mset_surf.vecs]
			rewrapz = array([[ionstraj[i,j,2]-1*(array(ionstraj[i,j,2])>mset_surf.vecs[j][2])*\
				mset_surf.vecs[j][2]+(array(ionstraj[i,j,2])<0.)*mset_surf.vecs[j][2] 
				for j in range(len(ionstraj[i]))] for i in range(len(ionstraj))])
			roundz = array([[int(rewrapz[i][j]/binws[j]) for j in range(len(rewrapz[i]))] 
				for i in range(len(rewrapz))])
			counts,bins = histogram(array(roundz).flatten(),range=(0,nzones),bins=nzones)
			#---Nb check the histogram with the following code
			edgebins = [[i.min()-1,i.max()+1] for i in [where((counts[:-1]/\
				float(sum(counts[:-1])))<ion_drop_thresh)[0]]][0]
			zonesub = range(edgebins[0]-sideslices+1,edgebins[0]+1)+range(edgebins[1],edgebins[1]+sideslices)
			print 'status: zonesub = '+str(zonesub)
			print 'status: requested bin width = '+str(bwid)
			print 'status: actual bin width with flush bins = '+str(mean(binws))
			print 'status: plotting histogram'
#			plt.grid(True);plt.hist(array(roundz).flatten(),range=(0,nzones),bins=nzones);plt.show()
			nzonessub = len(zonesub)
			#---Nb you could easily derive another zone selection and put it here

'''
ew = 80.10*8.854187 * 10**-12 #---F/m
kbt = 1.38*10**-23 * 300 #--m2kg/s2
zp = 1
zn = -1
np = len(ionsp[0])/product(mean(vecs,axis=0)[:3]) #---number per A3
nn = len(ionsn[0])/product(mean(vecs,axis=0)[:3]) #---number per A3
echarge = 1.6*10**-19
meter_to_ang = 10**10

lam = (4*pi/ew/kbt*echarge**2*(zp**2*np+zn**2*nn)*meter_to_ang)**2
lb = echarge**2/ew/kbt
lb = 8.1
def phi(r): 
	return kbt*lb/abs(r)*exp(-abs(r)/lam)
'''

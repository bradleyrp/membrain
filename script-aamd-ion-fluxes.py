#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = 'name P'

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
	'v509-40000-90000-10':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-ions',
		'trajsel':'s6-kraken-md.part0018.40000-90000-10.ions.xtc',
		'sysname_lookup_struct':'membrane-v509-atomP',
		'trajsel_struct':'s6-kraken-md.part0018.40000-90000-10.atomP.xtc',
		'ionname':'NA',
		'director':director_aamd_symmetric}}
analysis_names = [
	'v532-20000-58000-100',
	'v530-30000-100000-100',
	'v531-20000-62000-100',
	'v509-40000-90000-10'
	][-1:]
routine = ['compute',][1:0]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice
if 'compute' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		trajsel = trajsel_struct
		sysname_lookup = sysname_lookup_struct
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			keysysname='sysname_lookup_struct',keytrajsel='trajsel_struct')
		mset_surf = MembraneSet()
		mset_surf.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		#---collect ion positions
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		#---collect monolaye positions
		#---Nb this is where we can refine the method
		mset_surf.identify_monolayers(director)
		surf_select = mset_surf.universe.selectAtoms(selector)
		monoz = []
		for fr in whichframes:
			mset_surf.gotoframe(fr)
			surfpts = surf_select.coordinates()
			monoz.append([mean(surfpts[mset_surf.monolayer_residues[m],2],axis=0) for m in range(2)])
		#---construct bins		
		desired_binsize = 5
		upper_binws = [int(round((mset_surf.vec(i)[2]-monoz[i][0])/desired_binsize)) for i in whichframes]
		mid_binws = [int(round((monoz[i][0]-monoz[i][1])/desired_binsize)) for i in whichframes]
		lower_binws = [int(round((monoz[i][1])/desired_binsize)) for i in whichframes]
		#---old code frame-wise, mistaken because number of bins not consistent
		if 0:
			binedges = array([np.concatenate((
				linspace(0,monoz[i][0],lower_binws[i])[:-1],
				linspace(monoz[i][0],monoz[i][1],mid_binws[i])[:-1],
				linspace(monoz[i][1],mset_surf.vec(i)[2],upper_binws[i]))) 
				for i in whichframes])
			badframes = list(where([ionspos[j][:,2].max()>mset_surf.vec(j)[2] for j in range(5000)])[0])	
			whichframes2 = [i for i in whichframes if i not in badframes]
			disctraj  = array([[
				where(ionspos[j][:,2][i]<binedges[j][1:])[0][0] 
				for i in range(len(ionspos[0]))] for j in whichframes2]).T
		#---pick a standard number of bins for consistent zones
		bin_nums = [int(round(mean(i))) for i in [lower_binws,mid_binws,upper_binws]]
		badframes = list(where([ionspos[j][:,2].max()>mset_surf.vec(j)[2] for j in range(5000)])[0])	
		whichframes2 = [i for i in whichframes if i not in badframes]
		binedges = array([np.concatenate((
			linspace(0,monoz[i][1],bin_nums[0])[:-1],
			linspace(monoz[i][1],monoz[i][0],bin_nums[1])[:-1],
			linspace(monoz[i][0],mset_surf.vec(i)[2],bin_nums[2]))) 
			for i in whichframes])
		disctraj  = array([[
			where(ionspos[j][:,2][i]<binedges[j][1:])[0][0] 
			for i in range(len(ionspos[0]))] for j in whichframes2]).T
		#---quick calculation of a transition matrix
		disctrajt = disctraj.T
		tmat = zeros((len(binedges[0])-1,len(binedges[0])-1))
		for fr in range(len(disctrajt)-1):
			for i in range(len(disctrajt[fr])):
				tmat[disctrajt[fr][i],disctrajt[fr+1][i]] += 1
		#---plot the transition matrix
		plt.imshow(tmat,interpolation='nearest',origin='lower');plt.show()
		

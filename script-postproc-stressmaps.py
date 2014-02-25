#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

import scipy.interpolate
import scipy.integrate
from scipy import ndimage

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Key parameters
'''
Notes:

previously the parameters were set as follows
	framewise_test = [4,32,64,1]
	test = framewise_test
	span,nnum,numgridpts,distance_factor = test
	distance_factor was used in my old pseudo-RBF function
	nnum was the number of nearest neighbors that contributed to the sum smoothed out in the pseudo-RBF
	numgridpts is the 1D length of grid points sampled from the pseudo RBF
	
other than the span, which set the number of voxels integrated on both sides of the bilayer 
and the choice of how to define the midplane
the other parameters are all smoothing parameters
hence the way to calculate the c0 maps is to do the integration, dump the then-much-smaller data
and smooth it seperately
'''

#---analysis plan
analysis_descriptors = {
	'v614-120000-220000-200': 
		{'sysname':'membrane-v614',
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'nprots':4,
		'datdir3dpp':
			'/home/rpb/compbio/membrane-v614-enthx4-12800/a8-stress-s9-120000-220000/results',
		'framewise_part':4,
		'voxelsize':1.0},
	'v550-300000-400000-200': 
		{'sysname':'membrane-v550',
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'nprots':0,
		'datdir3dpp':
			'/home/rpb/compbio-alt/membrane-v550/'+\
			'a1-stress-1.0-framewise-md.part0006.300000-400000-200/results',
		'framewise_part':6,
		'voxelsize':1.0}}
analysis_names = ['v614-120000-220000-200','v550-300000-400000-200'][:]
routine = ['calc_c0maps']
span_sweep = [1,2,3,4,5,6]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'calc_c0maps' in routine:
	for aname in analysis_names:
		#---get names and load the structure
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		msetfile = 'pkl.structures.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		mset = unpickle(pickles+msetfile)
		picklename = 'pkl.stressmaps.'+specname_pickle(sysname,trajfile[0])+'pkl'
		result_data_spans = [MembraneData('collect_c0maps') for i in range(len(span_sweep))]
		for frame in range(len(mset.surf)):
			print 'status: computing curvature maps frame: '+str(frame)
			#---load raw stress tensor data
			file3dpp = pickles+'/'+datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+\
				str('%04d'%frame)+'.lp.dat3d'
			file3dpp = datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+str('%04d'%frame)+'.lp.dat3d'
			dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
			griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
			vecs = mset.vecs[frame]
			blocksize = (vecs/griddims)
			#---allocate space for the height-wise voxel tensions
			xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
				griddims[1]+2)])
			#---reference midplane surface
			ref_surf = mset.unzipgrid(mset.surf[frame],vecs=vecs)
			#---interpolate the reference surface to match the spacing of the stress tensors
			interp = scipy.interpolate.LinearNDInterpolator(ref_surf[:,0:2],ref_surf[:,2],fill_value=0.0)
			ref_surf_interp = array([[round((interp(i,j)+mset.surf_position[frame])/blocksize[2]) 
				for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
			#---loop over limits of integration
			for span in span_sweep:
				rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
				for pt in dat3dpp:
					if ((abs(pt[2]-ref_surf_interp[pt[0]][pt[1]]) <= span) and 
						(abs(pt[2]-ref_surf_interp[pt[0]][pt[1]]) <= span)):
						rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
							(pt[2]-ref_surf_interp[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
				#---integrate the collected stresses to determine the spontaneous curvature
				results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)]
					for i in range(griddims[0]+1)])
				#---set aside the data
				result_data_spans[span_sweep.index(span)].add(results,[frame])
		#---reorganize the data
		for r in range(len(result_data_spans)):
			result_data_spans[r].addnote(['c0_span',span_sweep[r]])
			for key in analysis_descriptors[aname]: 
				result_data_spans[r].addnote([key,analysis_descriptors[aname][key]])
			mset.store.append(result_data_spans[r])
		#---remove redundant data
		mset.xypts = []
		mset.protein = []
		mset.surf = []
		#---write the store, since saving mset is very costly for some reason
		#---Nb saving mset was 100MB with only 10 frames even after deleting the redundant data
		pickledump(mset.store,picklename,directory=pickles)

########### OLD JUNK JUNK JUNK

if mset.surf == [] and 0:
	aname = analysis_names[0]
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
	msetfile = 'pkl.structures.'+specname_pickle(sysname,trajfile[0])+'.pkl'
	mset = unpickle(pickles+msetfile)
	
if 0:
	vecs = mean(mset.vecs,axis=0)
	test = framewise_test
	logmaxminmean = []
	span,nnum,numgridpts,distance_factor = test
	allmaps = []
	#unzipsurfmean = mset.unzipgrid(list(mean(mset.surf,axis=0)),vecs=vecs)
	for frame in range(len(mset.surf))[::4]:
		print 'running frame = '+str(frame)
		#---Previously stored the frames in the pickle directory - now found with the simulations
		file3dpp = pickles+'/'+datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+\
			str('%04d'%frame)+'.lp.dat3d'
		file3dpp = datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+str('%04d'%frame)+'.lp.dat3d'
		dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
		griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
		vecs = mean(mset.vecs,axis=0)
		blocksize = (vecs/griddims)
		#---IMPORTANT CHANGE. I ADDED TRANSPOSE BELOW DUE TO MAJOR ERROR
		#---IMPORTANT CHANGE. NOW REMOVED AND FIXED UPSTREAM
		#unzipsurfmean = mset.unzipgrid([[mset.surf_position[frame] for j in range(mset.griddims[1])] for i in range(mset.griddims[0])],vecs=vecs)
		unzipsurfmean = mset.unzipgrid(mset.surf[frame],vecs=vecs)
		#unzipsurfmean = mset.unzipgrid(list(mean(mset.surf,axis=0)),vecs=vecs)
		#unzipsurfmean = mset.unzipgrid(list(array(mset.surf[frame]).T),vecs=vecs)
		themeanpos = mset.surf_position[frame]
		rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
		xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
			griddims[1]+2)])
		interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],
			fill_value=0.0)
		avginterpsurf = array([[round((interp(i,j)+themeanpos)/blocksize[2]) 
			for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
		for pt in dat3dpp:
			if ((abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and 
				(abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span)):
				rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
					(pt[2]-avginterpsurf[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
		results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)] 
			for i in range(griddims[0]+1)])
		allmaps.append(results)
if 0:
	#plt.imshow(array(results).T,interpolation='nearest',origin='lower');plt.show()
	plt.imshow((scipy.ndimage.filters.gaussian_filter(mean(allmaps,axis=0).T,3)),interpolation='nearest',origin='lower',cmap=mpl.cm.jet);plt.show()

	#smoothplot = mean([scipy.ndimage.filters.gaussian_filter(allmaps[i],1) for i in range(len(allmaps))],axis=0)
	#plt.imshow(smoothplot.T,interpolation='nearest',origin='lower',cmap=mpl.cm.jet);plt.show()

	
	
		

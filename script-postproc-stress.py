#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')

import os
import numpy as np
import scipy.interpolate
import scipy.integrate
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams.update({'font.style':'sans-serif'})
mpl.rcParams.update({'font.size': 16})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Parameters, initial (best)
nnnum = 5 			#---number of nearest neighbors
numgridpts = 15 	#---how many grid points, per direction, to sample (20 is optimal)
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
span = 3			#---how many voxels to move

#---Parameters, modified
voxelsize=1.0		#---voxel size from the rerun
flatref = False		#---use a flat reference midplane instead of the rounded average structure
view3d = False		#---view the map in 3D using mayavi
exag = 0 			#---exaggerates the peaks on the plot if you use the view3d
check_surf_mean = 0	#---whether to print the mean surface for 
protlenshow = 22 	#---Include only the first few protein points in the centroid calculation (o/w None)
plotunsmoothed = 0	#---Plot the raw data
brokeversion = 0	#---A weird effect with small spans and the wrong input to rezip grid
plotvisible = 0		#---Show the plots if you are doing the comparison
avgsurfround = 0	#---Whether you round the average surface to the voxel grid

#---Key parameters
framewise_test = [4,32,64,1]

#---Analysis plan
analysis_plan = slice(-1,None)
analysis_descriptors = [
	['v701.part0003.60000-160000-200','pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
		'localpressure.v701.part0003.60000-160000-200.3Dpp.dat',1,
		'/home/rpb/compbio-alt/membrane-v701-exo70-anti-dilute/a1-stress-1.0-framewise/results',
		'pkl.stressdecomp.membrane-v701.md.part0003.60000-160000-200.pkl',3],
	['v700.part0002.100000-200000-200','pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',
		'localpressure.v700.part0002.100000-200000-200.3Dpp.dat',1,
		'/home/rpb/compbio-alt/membrane-v700-exo70-dilute/a3-stress-1.0-framewise-100000-200000/results',
		'pkl.stressdecomp.membrane-v700.md.part0002.100000-200000-200.pkl',3],
	['v550.part0006.300000-400000-200','pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',
		'localpressure.v700.part0006.300000-400000-200.3Dpp.dat',0,
		'/home/rpb/compbio/membrane-v550/u5-stress-1.0-framewise-md.part0006.300000-400000-200/results',
		'pkl.stressdecomp.membrane-v550.md.part0006.300000-400000-200.pkl',6],
	['v614.part0002','pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl','',4,
		'/home/rpb/compbio/membrane-v614-enthx4-12800/a7-localpressure.v614.framewise/',
		'pkl.stressdecomp.membrane-v614-stress.md.part0002.rerun.pkl',2],
	['v612.part0003','pkl.structures.membrane-v612-stress.md.part0003.pkl','',2,
		'/home/rpb/compbio/membrane-v612-enthx1-12800/a3-localpressure.v612.framewise/',
		'pkl.stressdecomp.membrane-v612-stress.md.part0002.rerun.pkl',3],
	['v700.part0009.500000-700000-400','pkl.structures.membrane-v700.md.part0009.500000-700000-400.pkl',
        'localpressure.v701.part0009.500000-700000-400.3Dpp.dat',1,
        '/home/rpb/compbio-alt/membrane-v700-exo70-dilute/a4-stress-1.0-framewise-500000-700000-400/results',
        'pkl.stressdecomp.membrane-v700.md.part0009.500000-700000-400.pkl',9]]
		
#---Type of looping or parameter sweeps to do
#---Note: that I am temporarily dropping support for everything but batch_parameter_sweep_framewise
#---Note: that previous methods are in beta/script-postproc-stress-smooth-beta-backup-2014.01.08.py
batch_parameter_sweep_framewise = True

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def calculate_stressmap(span,nnnum,numgridpts,distance_factor,imagefile=None,runplots=True,logfile=None,
	plotunsmoothed=False,brokeversion=False):
	'''Curvature calculation via stress tensor, post-post processing.'''
	for pt in dat3dpp:
		if ((abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and 
			(abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span)):
			rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
				(pt[2]-avginterpsurf[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
	results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)] 
		for i in range(griddims[0]+1)])
	print shape(results)
	#---Pseudo-Radial basis function computation with K-nearest neighbors
	pts = mset.wrappbc(mset.unzipgrid(results/10**(-exag),vecs=mset.vec(0)),
		vecs=mset.vec(0),mode='grow',growsize=1.0)
	print shape(pts)
	#meshplot(pts,vecs=vecs,show='surf')
	#raw_input('...')
	newpts = array([[i,j] for i in linspace(-vecs[0],2*vecs[0],3*numgridpts) 
		for j in linspace(-vecs[1],2*vecs[1],3*numgridpts)])
	print shape(newpts)
	#---Distance matrix
	tree = scipy.spatial.cKDTree(pts[:,0:2])
	#---Unpack the points and perform the pseudo-RBF
	smoothed = []
	for pt in newpts:
		tmp = tree.query(pt,nnnum)
		#print tmp
		#print shape(tmp[0])
		smoothed.append([pt[0],pt[1],mean([pts[tmp[1][i],2]*1/((tmp[0][i])*sum(tmp[0]))**distance_factor 
			for i in range(len(tmp[0]))])])
	print 'shape(smoothed) = '+str(shape(smoothed))
	#nearset = [i for i in smoothed[1:] if np.linalg.norm(i[0:2]-mean(mset.protein[0],axis=0)[0:2])<300.]
	#---Drop infinite values and the edges, which are noisy
	nearset = [i for i in smoothed if np.isinf(i[2]) == False 
		and (i[0] != 0.0 and i[1] != 0.0 and i[0] != vecs[0] and i[1] != vecs[1])
		and (i[0] > 0. and i[0] < vecs[0] and i[1] > 0. and i[1] < vecs[1])]
	print 'shape(nearset) = '+str(shape(nearset))
	if view3d:
		meshplot(nearset,show='surf')
	if view3d and nprots > 0:
		meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,
			color=(0,0,0),opacity=0.2)
	#---Report maximum and mean curvatures assuming kappa = 20 kBT
	peakmax = (max(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	peakmin = (min(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	peakmean = (mean(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	print 'peakmax = '+str(peakmax)
	print 'peakmin = '+str(peakmin)
	print 'peakmean = '+str(peakmean)
	logmaxminmean.append([span,nnum,numgridpts,distance_factor,peakmax,peakmin,peakmean])
	global protlenshow
	if nprots > 0:
		protlen = int(shape(mset.protein[0])[0]/nprots)
		if protlenshow == None:
			protlenshow = protlen			
		protein_centers = mean([mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow] 
			for i in range(nprots)],axis=1)
	else:
		protein_centers = None
	if plotunsmoothed == True:
		dat = results
	elif brokeversion == True:
		dat = mset.rezipgrid(nearset,vecs=vecs)
	else:
		dat = mset.rezipgrid(nearset,vecs=vecs,grid=[numgridpts,numgridpts])
	return (dat,protein_centers)
	
def plot_stressmap(dat,protein_centers,nprots,numgridpts,imagefile=None,plotvisible=False):
	'''Plot a single stress map.'''
	fig = plt.figure()	
	ax = plt.subplot2grid((1,4),(0,0))
	plt.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft')
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=0.25,color='k',
				alpha=1.0)
			ax.add_patch(circ)
	if imagefile != None:
		plt.savefig(imagefile)
	if plotvisible:
		plt.show()
	else:
		plt.clf()

#---MAIN
#-------------------------------------------------------------------------------------------------------------
		
#---Perform the framewise analysis and dump to a datfile for reprocessing.
if batch_parameter_sweep_framewise:
	for ad in analysis_descriptors[analysis_plan]: 
		starttime = time.time()
		print 'Starting analysis job.'
		(systemname,msetfile,picklefile,nprots,datdir3dpp,framewise_out_name,framewise_part) = ad
		mset = unpickle(pickles+msetfile)
		test = framewise_test
		logmaxminmean = []
		span,nnum,numgridpts,distance_factor = test
		res_collection = []
		for frame in range(len(mset.surf)):
			print 'running frame = '+str(frame)
			#---Previously stored the frames in the pickle directory - now found with the simulations
			file3dpp = pickles+'/'+datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+str('%04d'%frame)+\
				'.lp.dat3d'
			file3dpp = datdir3dpp+'/'+'md.part'+str('%04d'%framewise_part)+'.fr'+str('%04d'%frame)+'.lp.dat3d'
			dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
			griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
			vecs = mean(mset.vecs,axis=0)
			blocksize = (vecs/griddims)
			#---IMPORTANT CHANGE. I ADDED TRANSPOSE BELOW DUE TO MAJOR ERRORRRRRR
			#unzipsurfmean = mset.unzipgrid(mset.surf[frame],vecs=vecs)
			unzipsurfmean = mset.unzipgrid(list(array(mset.surf[frame]).T),vecs=vecs)
			themeanpos = mset.surf_position[frame]
			rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
			xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
				griddims[1]+2)])
			interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],
				fill_value=0.0)
			avginterpsurf = array([[round((interp(i,j)+themeanpos)/blocksize[2]) 
				for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
			if ((themeanpos/blocksize[2] - span < 0.0) or 
				(themeanpos/blocksize[2] + span >  max(dat3dpp[:,2]))):
				print 'Warning: your span exceeds your box dimensions'
			res = calculate_stressmap(span,nnum,numgridpts,distance_factor,plotunsmoothed=False,
				brokeversion=False)
			res_collection.append(res)
			#plot_stressmap(res[0],res[1],nprots,numgridpts,imagefile=None,plotvisible=True)
		pickle.dump(res_collection,open(pickles+framewise_out_name,'w'))
		print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


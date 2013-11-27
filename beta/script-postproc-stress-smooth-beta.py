#!/usr/bin/python -i

newone = 1

if newone:
	from membrainrunner import *

	location = 'dark'
	execfile('locations.py')

	execfile('plotter.py')
	import scipy.interpolate
	import scipy.integrate


#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

'''
#---Parameters
nnnum = 5 			#---number of nearest neighbors
numgridpts = 20 	#---how many grid points, per direction, to sample (20 is optimal)
exag = -0 			#---exaggerates the peaks on the plot
distance_factor = 1	#---exponent on distance scaling in pseudo-RBF
check_surf_mean = 0	#---whether to print the mean surface for 
span = 5			#---How many voxels to move
voxelsize=1.0		#---Voxel size from the rerun
include_prots = 0	#---??????????

#---Parameters
nnnum = 5 			#---number of nearest neighbors
numgridpts = 40 	#---how many grid points, per direction, to sample (20 is optimal)
exag = -0 			#---exaggerates the peaks on the plot
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
check_surf_mean = 0	#---whether to print the mean surface for 
span = 10			#---How many voxels to move
voxelsize=1.0		#---Voxel size from the rerun
include_prots = 0	#---??????????

#---Parameters
nnnum = 5 			#---number of nearest neighbors
numgridpts = 15 	#---how many grid points, per direction, to sample (20 is optimal)
exag = -0 			#---exaggerates the peaks on the plot
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
check_surf_mean = 0	#---whether to print the mean surface for 
span = 3			#---How many voxels to move
voxelsize=1.0		#---Voxel size from the rerun
include_prots = 0	#---??????????

######### BEST BEST BEST
#---Parameters
nnnum = 5 			#---number of nearest neighbors
numgridpts = 15 	#---how many grid points, per direction, to sample (20 is optimal)
exag = -0 			#---exaggerates the peaks on the plot
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
check_surf_mean = 0	#---whether to print the mean surface for 
span = 3			#---How many voxels to move
voxelsize=1.0		#---Voxel size from the rerun
with_protein = 0	#---??????????

'''
#---Parameters
nnnum = 4 			#---number of nearest neighbors
numgridpts = 32 	#---how many grid points, per direction, to sample (20 is optimal)
exag = -1 			#---exaggerates the peaks on the plot
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
check_surf_mean = 0	#---whether to print the mean surface for 
span = 5			#---How many voxels to move
voxelsize=1.0		#---Voxel size from the rerun
with_protein = 0	#---??????????

'''
#---Batch settings
#---Order: nnum,numgridpts,distance_factor,span
param_batch = (
	[],
	[])
	
analysis_plan = (
	[],
	[])
'''

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if newone:
	mset = unpickle(pickles+'pkl.postproc-dimple-fit.membrane-v614.md.part0002.skip10.pkl')
	file3dpp = pickles+'localpressure.v614.part0002.3Dpp.dat'
	#mset = unpickle(pickles+'pkl.postproc-dimple-fit.membrane-v612.md.part0003.skip10.pkl')
	#file3dpp = pickles+'localpressure.v612.part0003.3Dpp.dat'
	#mset = unpickle(pickles+'pkl.structures.membrane-v550.md.part0008.rerun.pkl')
	#file3dpp = pickles+'localpressure.v550.part0008.3Dpp.dat'

#---Curvature calculation via stress tensor, post-post processing
#-------------------------------------------------------------------------------------------------------------

#---Original stress stuff ??????????
######## THIS IS WRONG, IT SUBTRACTS THE AVERAGE POSITION
dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
mset.calculate_average_surface()
vecs = mean(mset.vecs,axis=0)
#vecs=array([640,640,230])
griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],griddims[1]+2)])
unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],fill_value=0.0)
blocksize = (vecs/griddims)
themeanpos = mean(mset.surf_position)
#avginterpsurf = array([[round((interp(i,j)+mean(mset.surf_position))/blocksize[2])
#	for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
avginterpsurf = array([[themeanpos/blocksize[2]
	for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
if ((themeanpos/blocksize[2] - span < 0.0) or 
	(themeanpos/blocksize[2] + span >  max(dat3dpp[:,2]))):
	print 'Warning: your span exceeds your box dimensions'
for pt in dat3dpp:
	if ((abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and 
		(abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span)):
		rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
			(pt[2]-avginterpsurf[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)] 
	for i in range(griddims[0]+1)])

#---Pseudo-Radial basis function computation with K-nearest neighbors
pts = mset.wrappbc(mset.unzipgrid(results/10**exag,vecs=mset.vec(0)),
	vecs=mset.vec(0),mode='grow',growsize=1.5)		
newpts = array([[i,j] for i in linspace(0,vecs[0],numgridpts) for j in linspace(0,vecs[1],numgridpts)])

#---Distance matrix
tree = scipy.spatial.cKDTree(pts[:,0:2])

#---Unpack the points and perform the pseudo-RBF
smoothed = []
for pt in newpts:
	tmp = tree.query(pt,nnnum)
	smoothed.append([pt[0],pt[1],mean([pts[tmp[1][i],2]*1/((tmp[0][i])*sum(tmp[0]))**distance_factor 
		for i in range(len(tmp[0]))])])
#nearset = [i for i in smoothed[1:] if np.linalg.norm(i[0:2]-mean(mset.protein[0],axis=0)[0:2])<300.]

#---Drop infinite values and the edges, which are noisy
nearset = [i for i in smoothed if np.isinf(i[2]) == False 
	and (i[0] != 0.0 and i[1] != 0.0 and i[0] != vecs[0] and i[1] != vecs[1])]
meshplot(nearset,show='surf')
if with_protein:
	meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,color=(0,0,0),
		opacity=0.2)

#---Print maximum and mean curvatures assuming kappa = 20 kBT
print (max(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
print (min(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
print (mean(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)

if check_surf_mean:
	#---plot average structure
	meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,color=(0,0,0),
		opacity=0.5)
	meshplot(mset.surf_mean,vecs=mset.vec(0),show='surf')
	
#---


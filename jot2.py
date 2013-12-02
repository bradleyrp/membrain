#!/usr/bin/python

scan_span = [2,3,5,6]
scan_nnum = [4,16,32,64]
scan_numgridpts = [64,32,20,16]
scan_distance_factor = [2,1]

ordering = [scan_span,scan_nnum,scan_numgridpts,scan_distance_factor]
scans = [[i,j,k,l] for i in ordering[0] for j in ordering[1] for k in ordering[2] for l in ordering[3]]

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

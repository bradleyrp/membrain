#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')
import scipy.interpolate
import scipy.integrate

location = 'dark'
execfile('locations.py')

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

span = 10
voxelsize=1.0
tilesize=5

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---v614
if 1:
	mset = unpickle(pickles+'pkl.postproc-dimple-fit.membrane-v614.md.part0002.skip10.pkl')
	file3dpp = '/home/rpb/worker/worker-big/membrane-repository/membrane-v614-analyze/3Dpp-1.0.dat'
#---v032
if 0:
	mset = unpickle(pickles+'pkl.avgstruct.membrane-v032.md.part0002.skip10.pkl')
	file3dpp = '/home/rpb/worker/worker-big/membrane-repository/membrane-v032-exo70-antiparallel-revisit/3Dpp.dat'

dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
mset.calculate_average_surface()
vecs = mean(mset.vecs,axis=0)
griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],griddims[1]+2)])
unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],fill_value=0.0)
blocksize = (vecs/griddims)
avginterpsurf = array([[(interp(i,j)+mean(mset.surf_position))/blocksize[2] 
	for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
avginterpsurf = array([[mean(mset.surf_position)/blocksize[2]
	for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
if ((mean(mset.surf_position)/blocksize[2] - span < 0.0) or 
	(mean(mset.surf_position)/blocksize[2] + span >  max(dat3dpp[:,2]))):
	print 'Warning: your span exceeds your box dimensions'
for pt in dat3dpp:
	if ((abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and 
		(abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span)):
		rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
			(pt[2]-avginterpsurf[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)] 
	for i in range(griddims[0]+1)])

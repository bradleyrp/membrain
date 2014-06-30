#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')
import numpy

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

location = 'light'
execfile('locations.py')
#mset = unpickle(pickles+'pkl.avgstruct.membrane-v623-stress-test.md.part0005.skip10.pkl')
#mset = unpickle(pickles+'pkl.avgstruct.membrane-v614.md.part0002.skip10.pkl')
#mset = unpickle(pickles+'pkl.avgstruct.membrane-v599.relevant.pkl')
mset = unpickle(pickles+'pkl.avgstruct.membrane-v032.md.part0002.skip10.pkl')

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def curvcalc(z,lenscale):
	zy, zx  = numpy.gradient(z,lenscale)
	zxy, zxx = numpy.gradient(zx,lenscale)
	zyy, _ = numpy.gradient(zy,lenscale)
	H = (zx**2 + 1)*zyy - 2*zx*zy*zxy + (zy**2 + 1)*zxx
	H = -H/(2*(zx**2 + zy**2 + 1)**(1.5))
	K = ((zxx*zyy)-(zxy)**2)
	K = -K/(2*(zx**2 + zy**2 + 1)**(1.5))
	return [H,K]
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------
    
show_proteins = False
show_average = True
exaggeratem = 10**-1 #---reduce the heights to avoid shadows
exaggeratek = 10**0 #---reduce the heights to avoid shadows
    
#---Plot mean curvature
mset.calculate_average_surface()
protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]
if show_proteins:
	meshpoints(mset.protein[0]-[1*mean(mset.vecs,axis=0)[0],0,mean(mset.surf_position)-2.5],
		scale_factor=10,color=(0,0,0),opacity=0.5)
lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
cdat = lenscale*curvcalc(mset.surf[0],lenscale)[0]
curvs = []
for i in range(len(mset.surf)):
	curvs.append(curvcalc(mset.surf[i],lenscale)[0])
curvsm = mean(curvs,axis=0)
checkmesh(curvsm*exaggeratem,vecs=mean(mset.vecs,axis=0),tess=-2.5,wirecolor=(1,1,1),lsize=1,show='surf')
#---Plot gaussian curvature
if show_proteins:
	meshpoints(mset.protein[0]-[2*mean(mset.vecs,axis=0)[0],0,mean(mset.surf_position)-2.5],
		scale_factor=10,color=(0,0,0),opacity=0.5)
lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
curvsk = []
for i in range(len(mset.surf)):
	curvsk.append(curvcalc(mset.surf[i],lenscale)[1])
curvsmk = mean(curvsk*exaggeratek,axis=0)
checkmesh(curvsmk,vecs=mean(mset.vecs,axis=0),tess=-1.5,wirecolor=(1,1,1),lsize=1,show='surf')
print max([max(i) for i in curvsm])/200.
print min([min(i) for i in curvsm])/200.
if show_average:
	raw_input('save and close and then plot average...')
	if show_proteins:
		meshpoints(mset.protein[0]-[1*mean(mset.vecs,axis=0)[0],0,mean(mset.surf_position)-2.5],
			scale_factor=10,color=(0,0,0),opacity=0.5)
	checkmesh(mset.surf_mean,vecs=mean(mset.vecs,axis=0),tess=-1.5,wirecolor=(1,1,1),lsize=1,show='surf')
	

'''
oldcode

if 0:
	from membrainrunner import *
	execfile('plotter.py')
	mset = unpickle('pkl.avgstruct.membrane-v614.md.part0002.skip10.pkl')
	import numpy

if 0:
	from membrainrunner import *
	execfile('plotter.py')
	mset = unpickle('pkl.avgstruct.membrane-v700.md.part0005.skip10.half.pkl')
	import numpy

if 0:
	from membrainrunner import *
	execfile('plotter.py')
	import numpy
	pickles = '/store-delta/worker/worker-big/membrane-repository/pickle-repository/'
#	pickles = '/home/rpb/worker-big/membrane-repository/pickle-repository/'
#	mset = unpickle(pickles+'pkl.avgstruct.membrane-v612.md.part0003.skip10.pkl')
#	systemprefix = 'v612'
#	mset = unpickle(pickles+'pkl.avgstruct.membrane-v700.md.part0005.skip10.half.pkl')
#	systemprefix = 'v700'
#	mset = unpickle(pickles+'pkl.avgstruct.membrane-v623-stress-test.md.part0005.skip10.pkl')
	mset = unpickle(pickles+'pkl.postproc-dimple-fit.membrane-v614.md.part0002.skip10.protein2.pkl')
	mset.calculate_average_surface()

if 1:
	from membrainrunner import *
	execfile('plotter.py')
	import numpy
	pickles = '/home/rpb/worker-big/membrane-repository/pickle-repository/'
	#mset = unpickle(pickles+'pkl.avgstruct.membrane-v623-stress-test.md.part0005.skip10.pkl')
	#mset = unpickle(pickles+'pkl.avgstruct.membrane-v614.md.part0002.skip10.pkl')
	mset = unpickle(pickles+'pkl.avgstruct.membrane-v599.relevant.pkl')
	mset.calculate_average_surface()

def curvcalc(z,lenscale):
	zy, zx  = numpy.gradient(z,lenscale)
	zxy, zxx = numpy.gradient(zx,lenscale)
	zyy, _ = numpy.gradient(zy,lenscale)
	H = (zx**2 + 1)*zyy - 2*zx*zy*zxy + (zy**2 + 1)*zxx
	H = -H/(2*(zx**2 + zy**2 + 1)**(1.5))
	K = ((zxx*zyy)-(zxy)**2)
	K = -K/(2*(zx**2 + zy**2 + 1)**(1.5))
	return [H,K]
    
if 0:
	Z = mset.surf[0]
	X = mset.rezipgrid(mset.unzipgrid(mset.surf[0],vecs=mset.vecs[0]),whichind=0)
	Y = mset.rezipgrid(mset.unzipgrid(mset.surf[0],vecs=mset.vecs[0]),whichind=1)
	xyztmp = transpose([X,Y,Z])
	xyz = swapaxes(xyztmp,1,2)
	cdat2 = curvcalc(xyz)[0]
    
if 1:
	#protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]
	#meshpoints(mset.protein[0]-[0,0,mean(mset.surf_position)-2.5],scale_factor=20,color=(1,1,1))
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	cdat = lenscale*curvcalc(mset.surf[0],lenscale)[0]
	curvs = []
	for i in range(len(mset.surf)):
		curvs.append(200*curvcalc(mset.surf[i],lenscale)[0])
	curvsm = mean(curvs,axis=0)
	#meshplot(curvsm,vecs=mean(mset.vecs,axis=0))
	checkmesh(curvsm,vecs=mean(mset.vecs,axis=0),tess=-1.5,wirecolor=(1,1,1),lsize=1,show='surf')
	#raw_input()
	print max([max(i) for i in curvsm])/200.
	print min([min(i) for i in curvsm])/200.
	
if 1:
	#protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]
	#meshpoints(mset.protein[0]-[0,0,mean(mset.surf_position)-2.5],scale_factor=20,color=(1,1,1))
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	curvsk = []
	for i in range(len(mset.surf)):
		curvsk.append(curvcalc(mset.surf[i],lenscale)[1])
	curvsmk = mean(curvsk,axis=0)
	checkmesh(curvsmk,vecs=mean(mset.vecs,axis=0),tess=-2.5,wirecolor=(1,1,1),lsize=1,show='surf')
	#meshplot(curvsmk,vecs=mean(mset.vecs,axis=0))

if 0:
	mset.calculate_average_surface()
	meshplot(mset.surf_mean,vecs=mean(mset.vecs,axis=0),opacity=1.,show='both')
	protcom = mean(mset.protein[0],axis=0)-[0,0,mean(mset.surf_position)]-2.5
	meshpoints(mset.protein[0]-[0,0,mean(mset.surf_position)-2.5],scale_factor=20,color=(1,1,1))
'''

#!/usr/bin/python

from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq

def g2dresid(params,x,y,z):
	a,b,c1,c2,w1,w2,t1 = params
	return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)-z

def g2d(params,x,y):
	a,b,c1,c2,w1,w2,t1 = params
	return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)

if 0:
	#fr = 10
	#mset.gotoframe(fr)
	prot_nres = 158
	#mset.calculate_midplane(selector,fr)
	protpts = mset.universe.selectAtoms('name BB').coordinates()
	meshplot(mset.surf[2]+mset.surf_position[2],vecs=mset.vec(fr))
	meshpoints(protpts,scale_factor=10)
if 0:
	mset.calculate_average_surface()
	fr = mset.surf_index[5]
	mset.gotoframe(fr)
	protpts = mset.universe.selectAtoms('name BB').coordinates()
	meshpoints(mean(protpts,axis=0),scale_factor=50)
	#meshplot(mset.surf[mset.surf_index.index(fr)]+mset.surf_position[mset.surf_index.index(fr)],vecs=mset.vec(fr))
	#meshplot(mset.surf_mean+mean(mset.surf_position,axis=0),vecs=mset.vec(fr))
	meshplot(mset.surf_mean,vecs=mset.vecs[0],maxmin=[-10.,10.])
if 0:
	mean_params = mean(gauss_params,axis=0)
	mean_params = gauss_params[12]
	surfpts = array([[i[0],i[1],g2d(mean_params,i[0],i[1])-mean(top,axis=0)[2]] for i in near])
	#meshplot(surfpts,show='surf')
	#X,Y,Z = surfpts[:,0],surfpts[:,1],surfpts[:,2]
	#X,Y = meshgrid(linspace(center[0]-100,center[0]+100,20),linspace(center[1]-100,center[1]+100,20))
	#X = linspace(0,mset.vecs[0][0],10)
	#Y = linspace(0,mset.vecs[0][1],10)
	xypts = array([[i,j] for i in linspace(0,mset.vecs[0][0],10) for j in linspace(0,mset.vecs[0][1],10)])
	X = xypts[:,0]
	Y = xypts[:,1]
	Z = array([[g2d(gauss_params[0],X[i],Y[j]) for i in range(len(X))] for j in range(len(Y))] )
	mysurf = mlab.contour_surf(X,Y,Z)
	raw_input()
	
if 0:
	for j in range(100):
		surfpts = array([[i[0],i[1],g2d(gauss_params[j%len(gauss_params)],i[0],i[1])-mean(top,axis=0)[2]] for i in near])
		#X,Y,Z = surfpts[:,0],surfpts[:,1],surfpts[:,2]
		Z = array([g2d(gauss_params[j%len(gauss_params)],X[i],Y[i]) for i in range(len(X))])
		#X,Y=meshgrid(linspace(center[0]-10,center[0]+10,10),linspace(center[1]-10,center[1]+10,10))
		#mysurf.mlab_source.x = X
		#mysurf.mlab_source.y = Y
		mysurf.mlab_source.z = Z
		raw_input()
if 0:
	gauss_params = []
	for fr in range(len(mset.universe.trajectory)):
		mset.gotoframe(fr)
		top = array([mean(mset.universe.residues[i].selectAtoms('name PO4').coordinates(),axis=0) for i in mset.monolayer_residues[0]])
		protpts = mset.universe.selectAtoms('name BB').coordinates()
		center = mean(protpts,axis=0)
		near = array([i for i in top if norm(i-center) < 50])
		#meshplot(near-[0,0,mean(top,axis=0)[2]],vecs=mset.vec(fr),maxmin=[-10.,10.])
		#meshpoints(center-[0,0,mean(top,axis=0)[2]],scale_factor=50)
		p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),args=(near[:,0],near[:,1],near[:,2]))
		print 'frame = '+str(fr)
		print list(p_opt[0])
		gauss_params.append(p_opt[0])
		surfpts = array([[i[0],i[1],g2d(p_opt[0],i[0],i[1])-mean(top,axis=0)[2]] for i in near])
		#meshplot(surfpts,show='surf')
	gauss_params = array(gauss_params)
if 1:
	mysurf = mlab.surf(array([[g2d(gauss_params[1],x,y) for x in linspace(0,mset.vecs[0][0],10)] for y in linspace(0,mset.vecs[0][1],10)]))
	for k in range(100):
		Z = array([[g2d(gauss_params[k%len(gauss_params)],x,y) for x in linspace(0,mset.vecs[0][0],10)] for y in linspace(0,mset.vecs[0][1],10)])
		mysurf.mlab_source.scalars = Z
		print Z[0][0]
		raw_input()
		

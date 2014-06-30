#!/usr/bin/python -i 

#from membrainrunner import *
#execfile('plotter.py')
#mset = unpickle('pkl.postproc-dimple.v612.pkl')

from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq


def g2dresid(params,x,y,z):
	a,b,c1,c2,w1,w2,t1 = params
	return 0+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)-z

def g2d(params,x,y):
	a,b,c1,c2,w1,w2,t1 = params
	return 0+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)

if 0:
	meshplot(mset.surf_mean,vecs=mset.vec(0))
	meshpoints(mean(mset.protein,axis=0)-[0,0,mean(mset.surf_position)],scale_factor=50)
	
if 0:
	meshpoints(mean(mset.protein,axis=0)-[0,0,mean(mset.surf_position)],scale_factor=20,color=(1,1,1))
	meshplot(mset.surf_mean,vecs=mset.vec[0])
	gauss_params = []
	#for fr in mset.surf_index[0:10]:
	for fr in [0]:
		print 'Fitting Gaussian for frame '+str(fr)
		pts = mset.unzipgrid(mset.surf[mset.surf_index.index(fr)],vecs=mset.vec(fr))
		protpts = mset.protein[mset.surf_index.index(fr)][0:158]
		#meshplot(pts,vecs=mset.vec(fr),show='wire')
		#meshpoints(mean(mset.protein,axis=0)-[0,0,mean(mset.surf_position)],scale_factor=20,color=(1,1,1))
		protcom = mean(protpts,axis=0)-[0,0,mean(mset.surf_position)]
		#meshpoints(protcom,scale_factor=50,color=(0,0,0))
		near = array([i for i in pts if norm(i[0:2]-protcom[0:2]) < 100])
		#meshplot(near)
		center = protcom
		p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),args=(near[:,0],near[:,1],near[:,2]))
		gauss_params.append(p_opt[0])
		guess = array([[i[0],i[1],g2d(gauss_params[-1],i[0],i[1])] for i in pts if norm(i[0:2]-protcom[0:2]) < 200])
		meshplot(guess,opacity=0.2)

# retry

#--original
if 0:
	gauss_params = []
	domainw = 2
	fr = 4
	protpts = mset.protein[mset.surf_index.index(fr)][0:158]
	protcom = mean(protpts,axis=0)-[0,0,mean(mset.surf_position)]
	domain,domain2 = [coarse_data[fr][1],coarse_data[fr][2]]
	meshplot(mset.surf[0])
	#meshpoints(mean(mset.protein,axis=0)-[0,0,mean(mset.surf_position)-30],scale_factor=20,color=(1,1,1))
	points = array([[i[0],i[1],i[2]] for i in  coarse_data[fr][0]])
	ptsfilt = array([i for i in coarse_data[fr][0] if (domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 and domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1)])
	ptsfilt = array([[i[0],i[1],i[2]] for i in ptsfilt])
	meshplot(array([[i[0],i[1],10*i[2]] for i in ptsfilt]))
	center = [mean(ptsfilt[:,0]),mean(ptsfilt[:,1])]
	p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),args=(ptsfilt[:,0],ptsfilt[:,1],ptsfilt[:,2]))
	gauss_params.append(p_opt[0])
	guess = array([[i[0],i[1],10*g2d(gauss_params[-1],i[0],i[1])] for i in points[:,0:2]])
	meshplot(guess,opacity=0.4)

#--trying to scale right
if 0:
	gauss_params = []
	domainw = 2
	fr = 4
	protpts = mset.protein[mset.surf_index.index(fr)][0:158]
	protcom = mean(protpts,axis=0)-[0,0,mean(mset.surf_position)]
	domain,domain2 = [coarse_data[fr][1],coarse_data[fr][2]]
	lenscale = (mset.vecs[0]/(mset.griddims[0]+1))[0]
	meshpoints(mean(mset.protein,axis=0)-[0,0,mean(mset.surf_position)-30],scale_factor=20,color=(1,1,1))
	points = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in  coarse_data[fr][0]])
	ptsfilt = array([i for i in coarse_data[fr][0] if (domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 and domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1)])
	ptsfilt = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in ptsfilt])
	meshplot(array([[i[0],i[1],10*i[2]] for i in ptsfilt]))
	center = [mean(ptsfilt[:,0]),mean(ptsfilt[:,1])]
	p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),args=(ptsfilt[:,0],ptsfilt[:,1],ptsfilt[:,2]))
	gauss_params.append(p_opt[0])
	guess = array([[i[0],i[1],10*g2d(gauss_params[-1],i[0],i[1])] for i in points[:,0:2]])
	meshplot(guess,opacity=0.4)

# current
# get camera position
if 1:
	#cam = mlab.move()
	fr = 1
	meshpoints(mset.protein[fr]-[0,0,mean(mset.surf_position)+25.],scale_factor=20,color=(1,1,1))
	meshplot(mset.surf[0],vecs=mset.vecs[0])
	v=mlab.view()
	boxcenter = (mean(mset.vecs,axis=0)/2.)[0:2]
	domainw = 2
	gauss_params = []
	raw_input("Select the desired camera angle.")
if 1:
	v=mlab.view()
	mlab.clf()
	#mlab.options.offscreen = True
	for fr in range(10):
		print 'processing frame '+str(fr)
		domain,domain2 = [coarse_data[fr][1],coarse_data[fr][2]]
		meshpoints(mset.protein[fr]-[0,0,mean(mset.surf_position)+25.],scale_factor=20,color=(1,1,1))
		lenscale = (mset.vecs[0]/(mset.griddims[0]+1))[0]
		points = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in  coarse_data[fr][0]])
		ptsfilt = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in coarse_data[fr][0] if (domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 and domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1)])
		#meshplot(points,show='wire')
		#meshplot(ptsfilt)
		meshpoints(ptsfilt,scale_factor=10)
		center = [mean(ptsfilt[:,0]),mean(ptsfilt[:,1])]
		p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),args=(ptsfilt[:,0],ptsfilt[:,1],ptsfilt[:,2]))
		gauss_params.append(p_opt[0])
		guess = array([[i[0],i[1],g2d(gauss_params[-1],i[0],i[1])] for i in points[:,0:2] if norm(i[0:2]-boxcenter)<200])
		meshplot(guess,opacity=0.5,show='surf')
		meshplot(mset.surf[fr],vecs=mset.vec(fr),show='wire',opacity=0.2,wirecolor=(1,1,1))
		#mlab.move(cam[0])
		#mlab.view(focalpoint=cam[1])
		mlab.view(*v)
		filename=('mlab-fig%05d.png'%fr)
		mlab.savefig(filename,size=(1920,1080))
		mlab.clf()
		mlab.close()
	
	
	
	
	
	

############################################

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
if 0:
	mysurf = mlab.surf(array([[g2d(gauss_params[1],x,y) for x in linspace(0,mset.vecs[0][0],10)] for y in linspace(0,mset.vecs[0][1],10)]))
	for k in range(100):
		Z = array([[g2d(gauss_params[k%len(gauss_params)],x,y) for x in linspace(0,mset.vecs[0][0],10)] for y in linspace(0,mset.vecs[0][1],10)])
		mysurf.mlab_source.scalars = Z
		print Z[0][0]
		raw_input()
		

#!/usr/bin/python -i 

#from membrainrunner import *
#execfile('plotter.py')
#mset = unpickle('pkl.postproc-dimple.v612.pkl')

from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq


def g2dresid(params,x,y,z):
	a,b,c1,c2,w1,w2,t1 = params
	return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)-z

def g2d(params,x,y):
	a,b,c1,c2,w1,w2,t1 = params
	return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)

######## why is everything here in angstroms?

if 1:
	domainw = 2
	gauss_params = []
	boxcenter = (mean(mset.vecs,axis=0)/2.)[0:2]
	fr = 0
	meshpoints(mset.protein[fr]-[0,0,mean(mset.surf_position)+25.],scale_factor=20,color=(1,1,1))
	meshplot(mset.surf[0],vecs=mset.vecs[0],show='wire')
	raw_input("Select the desired camera angle.")
	v=mlab.view()
	mlab.clf()
	for fr in range(len(coarse_data)):
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
	
# print things for mathematica
if 1:
	for fr in range(len(coarse_data)):
		print 'Writing frame to file '+str(fr)
		surfzunfold = coarse_data[fr][0]
		domain = coarse_data[fr][1]
		protptsplus = coarse_data[fr][2]
		domainw = coarse_data[fr][3]
		frameno = coarse_data[fr][4]
		write_framefilter(surfzunfold,array(domain),array(protptsplus),domainw,frameno)

if 1:
	filept = open(systemprefix+'-framefilter/params.dat','w')
	for i in gauss_params:
		for j in i:
			filept.write(str(j)+'\t')
		filept.write('\n')
	filept.close()





		

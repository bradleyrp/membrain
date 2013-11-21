#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')
import numpy as N
import pylab
from scipy.optimize import curve_fit
from numpy import array
from scipy.optimize import leastsq
import os

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
skip = 1
framecount = None
location = 'light'
execfile('locations-rpb.py')

#---Settings
do_framefilter =0
write_framefilters = 0
dofits = 0
write_params = 0
make_video = 0
#---filenames for different runs
systemprefix_in = 'membrane-v614.md.part0002.skip10'
#systemprefix_in = 'membrane-v599.relevant'
#systemprefix = systemprefix_in+'.protein-bigneighbor'
systemprefix = systemprefix_in+'.standard'
#---end filenames for different runs
startpickle = pickles+'pkl.avgstruct.'+systemprefix_in+'.pkl'
#---Slice of the protein points to use in the filter
protein_subset_slice = slice(None)
#---Specify which direction (1,-1) you are seeking curvature
heightfilter = -1
#---Tile/filter size
domainw = 4
#---How many manhattan lengths around the protein to include in the filter
manhatdist = 1

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def curvature_coarsen(slc,prot,domainw,frameno,rounder=20.,write=False):
	'''Take a segment of the bilayer near the protein and save it for fitting.'''
	prot = [[i[0]/rounder,i[1]/rounder,i[2]/rounder] for i in prot]
	m = shape(slc)[0]
	n = shape(slc)[1]
	#---Unfold
	surfzunfold = []
	for i in range(shape(slc)[0]):
		for j in range(shape(slc)[1]):
			surfzunfold.append([i,j,slc[i,j]])
	#---Coarsen 
	compacted = []
	for i in range(len(surfzunfold)):
		compacted.append([int(round(j/domainw)) for j in surfzunfold[i][0:2]]+[surfzunfold[i][2]])
	compacted = array(compacted)
	#---Group blocks
	mc = int(max(compacted[:,0])+1)
	nc = int(max(compacted[:,1])+1)
	averaged = [[[] for i in range(nc)] for j in range(mc)]
	for i in compacted:
		averaged[int(i[0])][int(i[1])].append(i[2])
	#---Average the groups
	for i in range(mc):
		for j in range(nc):
			averaged[i][j] = mean(averaged[i][j])
	#---Filter heights
	for i in range(mc):
		for j in range(nc):
			if (averaged[i][j] < 0 and heightfilter == 1) or (averaged[i][j] > 0 and heightfilter == -1):
				averaged[i][j] = 1
			else:
				averaged[i][j] = 0
	#---Scale proteins
	protcopy = list(array(list(prot)))
	for i in range(shape(protcopy)[0]):
		for j in range(2):
			protcopy[i][j] = int(round(protcopy[i][j]/domainw))
	protpts = [[0 for i in range(nc)] for j in range(mc)]
	protptsplus = [[0 for i in range(nc)] for j in range(mc)]
	#---Redundant frame border
	for i in protcopy:
		protpts[int(i[0])][int(i[1])] = 1
		for xd in range(-manhatdist,manhatdist+1):
			for yd in range(-manhatdist,manhatdist+1):
				protptsplus[int(i[0])+xd][int(i[1])+yd] = 1
	cenx1 = int(round(shape(averaged)[0]*3/4))
	ceny1 = int(round(shape(averaged)[0]*1/4))
	cenx2 = int(round(shape(averaged)[0]*1/4))
	ceny2 = int(round(shape(averaged)[0]*3/4)-2)
	domain = [[0 for i in range(nc)] for j in range(mc)]
	domain[cenx1][ceny1] = 1
	domain[cenx2][ceny2] = 1
	neighbors = []
	for i in protcopy:
		if protpts[int(i[0])][int(i[1])] == 1:
			neighbors.append(i)
		while len(neighbors) > 0:
			pt = [int(k) for k in neighbors.pop(0)]
			if averaged[pt[0]][pt[1]] == 0:
				domain[pt[0]][pt[1]] = 1
				if pt[0] < mc-1 and domain[pt[0]+1][pt[1]] == 0:
					neighbors.append([pt[0]+1,pt[1]])
				if pt[0] > 0 and domain[pt[0]-1][pt[1]] == 0:
					neighbors.append([pt[0]-1,pt[1]])
				if pt[1] < nc-1 and domain[pt[0]][pt[1]+1] == 0:
					neighbors.append([pt[0],pt[1]+1])
				if pt[1] > 0 and domain[pt[0]][pt[1]-1] == 0:
					neighbors.append([pt[0],pt[1]-1])
			else:
				domain[pt[0]][pt[1]] = -1
	if write:
		write_framefilter(surfzunfold,array(domain),array(protptsplus),domainw,frameno)
	hinton_custom(array(averaged),array(protpts),array(protptsplus),array(domain),
		filename=(pickles+systemprefix+'-framefilter/figs'+'/fig%05d.png'%frameno))
	return [surfzunfold,array(domain),array(protptsplus),domainw,frameno]
		
def write_framefilter(slc,domain,domain2,domainw,frameno,lenscale):
	'''Write the filtered/coarsened frame for use in another plotting program, namely Mathematica.'''
	filpt = open(pickles+systemprefix+'-framefilter/framefilter.%05d.dat'%frameno,'w')
	for i in slc:
		if domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 and \
			domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1:
			filpt.write(str(lenscale*i[0])+' '+str(lenscale*i[1])+' '+str(i[2])+'\n')
	filpt.close()

def batch_framefilter(dwsize):
	'''Batch caculate and write filtered/coarsened frames.'''
	result_data = MembraneData('dimple_filter')
	for i in range(0,len(mset.surf),10): ########### hacked
		print 'Computing framefilter for frame: '+str(i)
		result_data.add(curvature_coarsen(mset.surf[i],proteins[i],dwsize,i),[i])
	mset.store.append(result_data)

def view_average():
	'''Plot the average structure.'''
	mset.calculate_average_surface()
	meshplot(mset.surf_mean,vecs=mset.vecs[0])
	protpts = array([i-[0,0,25] for i in mean(mset.protein,axis=0)])
	meshpoints(protpts-[0,0,mean(mset.surf_position)],scale_factor=10,color=(1,1,1))
	
def g2d(params,x,y):
	'''Fitting function.'''
	a,b,c1,c2,w1,w2,t1 = params
	t1 = pi/4
	a = 0
	return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	
def g2dresid(params,x,y,z):
	'''Fitting function residual.'''
	return g2d(params,x,y)-z

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Make new directories for old-fashioned output and figure storage
if not os.path.isdir(pickles+systemprefix+'-framefilter'):
	os.mkdir(pickles+systemprefix+'-framefilter')
else:
	print 'Beware possible overwrites in '+pickles+systemprefix+'-framefilter'
if not os.path.isdir(pickles+systemprefix+'-framefilter/figs'):
	os.mkdir(pickles+systemprefix+'-framefilter/figs')
else:
	print 'Beware possible overwrites in '+pickles+systemprefix+'-framefilter/figs'
if not os.path.isdir(pickles+systemprefix+'-framefilter/figs-mlab'):
	os.mkdir(pickles+systemprefix+'-framefilter/figs-mlab')
else:
	print 'Beware possible overwrites in '+pickles+systemprefix+'-framefilter/figs-mlab'

#---Filter the frames if necessary
if do_framefilter:
	mset = unpickle(startpickle)
	proteins_all = array(mset.protein)
	proteins = proteins_all[:,protein_subset_slice]
	print 'Filtering frames.'
	batch_framefilter(domainw)
	print 'Saving pickle'
	pickledump(mset,pickles+'pkl.postproc-dimple-fit.'+systemprefix+'.pkl')
else:
	mset = unpickle(pickles+'pkl.postproc-dimple-fit.'+systemprefix+'.pkl')

#---Write the frame-filtered data
if write_framefilters:
	data = mset.getdata('dimple_filter').data
	for fr in range(len(data)):
		print 'Writing frame to file '+str(fr)
		surfzunfold = data[fr][0]
		domain = data[fr][1]
		protptsplus = data[fr][2]
		domainw = data[fr][3]
		frameno = data[fr][4]
		lenscale = (mset.vec(fr)/(mset.griddims[0]))[0]
		write_framefilter(surfzunfold,array(domain),array(protptsplus),domainw,frameno,lenscale)

#---Fitting procedure
if dofits:
	print 'Starting the fitting procedure.'
	gauss_params = []
	boxcenter = (mean(mset.vecs,axis=0)/2.)[0:2]
	for fr in range(0,len(mset.surf),20): #########33 hacked
		print 'processing frame '+str(fr)
		points_original = mset.getdata('dimple_filter').get(['frame',fr,'type','points'])	
		domain = mset.getdata('dimple_filter').get(['frame',fr,'type','domain'])
		domain2 = mset.getdata('dimple_filter').get(['frame',fr,'type','protptsplus'])
		lenscale = (mset.vecs[0]/(mset.griddims[0]+1))[0]
		points = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in  points_original])
		ptsfilt = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in points_original 
			if (domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 
			and domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1)])
		if len(ptsfilt) == 0:
			print 'No filtered points here'
			gauss_params.append([0,1,center[0],center[1],50,50,0])
		else:
			center = [mean(ptsfilt[:,0]),mean(ptsfilt[:,1])]
			p_opt = leastsq(g2dresid,array([0,1,center[0],center[1],50,50,0]),
				args=(ptsfilt[:,0],ptsfilt[:,1],ptsfilt[:,2]))
			gauss_params.append(p_opt[0])
		

#---Old-fashioned saving for mathematica
if write_params:
	filept = open(pickles+systemprefix+'-framefilter/params.dat','w')
	for i in gauss_params:
		for j in i:
			filept.write(str(j)+'\t')
	filept.close()
	
#---Make video
if make_video:
	gauss_params = [[float(i) for i in line.strip().split()] 
		for line in open(pickles+systemprefix+'-framefilter/params.dat')]
	fr = 0
	meshpoints(mset.protein[fr]-[0,0,mean(mset.surf_position)+25.],scale_factor=20,color=(1,1,1))
	meshplot(mset.surf[0],vecs=mset.vecs[0],show='wire')
	raw_input("Select the desired camera angle.")
	v=mlab.view()
	mlab.clf()
	boxcenter = (mean(mset.vecs,axis=0)/2.)[0:2]
	for fr in range(len(mset.surf)):
		print 'Printing frame '+str(fr)
		meshpoints(mset.protein[fr]-[0,0,mean(mset.surf_position)+25.],scale_factor=20,color=(1,1,1))
		points_original = mset.getdata('dimple_filter').get(['frame',fr,'type','points'])	
		domain = mset.getdata('dimple_filter').get(['frame',fr,'type','domain'])
		domain2 = mset.getdata('dimple_filter').get(['frame',fr,'type','protptsplus'])
		lenscale = (mset.vecs[0]/(mset.griddims[0]+1))[0]
		points = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in  points_original])
		ptsfilt = array([[i[0]*lenscale,i[1]*lenscale,i[2]] for i in points_original 
			if (domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 
			and domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1)])
		meshpoints(ptsfilt,scale_factor=10)
		center = [mean(ptsfilt[:,0]),mean(ptsfilt[:,1])]
		guess = array([[i[0],i[1],g2d(gauss_params[fr],i[0],i[1])] for i in points[:,0:2] 
			if norm(i[0:2]-boxcenter)<200])
		meshplot(guess,opacity=0.5,show='surf')
		meshplot(mset.surf[fr],vecs=mset.vec(fr),show='wire',opacity=0.2,wirecolor=(1,1,1))
		mlab.view(*v)
		filename=('mlab-fig%05d.png'%fr)
		mlab.savefig(pickles+systemprefix+'-framefilter/figs-mlab/'+filename,size=(1920,1080))
		mlab.clf()
		mlab.close()


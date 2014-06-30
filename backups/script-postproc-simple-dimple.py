#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')
import numpy as N
import pylab

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

manhatdist = 1
mset = unpickle('pkl.avgstruct.membrane-v614.md.part0002.skip10.pkl')
systemprefix = 'v614'
#mset = unpickle('pkl.avgstruct.membrane-v612.md.part0003.skip10.pkl')
#systemprefix = 'v612'
surfaces = mset.surf
proteins = array([i[0:158] for i in mset.protein])
heightfilter = 1 # choose which direction to filter

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def hinton_blob(x,y,area,colour):
    hs = numpy.sqrt(area) / 2                                                                                           
    xcorners = numpy.array([x - hs, x + hs, x + hs, x - hs])                                                            
    ycorners = numpy.array([y - hs, y - hs, y + hs, y + hs])                                                            
    pylab.fill(xcorners, ycorners, colour, edgecolor=colour,alpha=0.5)                                                            
                                                                                                                    
def hinton_custom(W,protblob,protblobplus,domain,maxWeight=None,show=False,filename=None):
	reenable = False                                                                                               
	if pylab.isinteractive():                                                                                           
		pylab.ioff()
		reenable = True                                                                                                 
	fig1 = plt.gcf()
	pylab.clf()                                                                                                         
	height, width = W.shape                                                                                         
	if not maxWeight:                                                                                               
		maxWeight = 2**numpy.ceil(numpy.log(numpy.max(numpy.abs(W)))/numpy.log(2))                                                      
	pylab.fill(numpy.array([0,width,width,0]),numpy.array([0,0,height,height]),'gray')                                          
	pylab.axis('off')                                                                                                 
	pylab.axis('equal')                                                                                                 
	for x in xrange(width):                                                                                         
		for y in xrange(height):
			_x = x+1
			_y = y+1
			w = W[y,x]
			if w == 0:
				color = '#FFFF33'
			elif w == 1:
				color = '#0066CC' #lightish   
			hinton_blob(_x - 0.5, height - _y + 0.5, 1,color)
			if domain[y,x] == 1:
				hinton_blob(_x - 0.5, height - _y + 0.5, 1,'#FF9900')
			if protblobplus[x,y] == 1:
				hinton_blob(_x - 0.5, height - _y + 0.5, 0.4,'WHITE')
			if protblob[x,y] == 1:
				hinton_blob(_x - 0.5, height - _y + 0.5, 0.25,'BLACK')
	if reenable:
		pylab.ion()
	if filename != None:
		fig1.savefig(filename,dpi=200)
	if show:
		pylab.show()
	else:
		pylab.clf()

def write_framefilter(slc,domain,domain2,domainw,frameno):
	filpt = open(systemprefix+'-framefilter/framefilter.%05d.dat'%frameno,'w')
	for i in slc:
		if domain[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1 and \
			domain2[int(round(i[0]/domainw))][int(round(i[1]/domainw))] == 1:
			filpt.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
	filpt.close()
	
def coarsen3(slc,prot,domainw,frameno,rounder=20.):
	prot = [[i[0]/rounder,i[1]/rounder,i[2]/rounder] for i in prot]
	# get shape
	m = shape(slc)[0]
	n = shape(slc)[1]
	# unfold
	surfzunfold = []
	for i in range(shape(slc)[0]):
		for j in range(shape(slc)[1]):
			surfzunfold.append([i,j,slc[i,j]])
	# coarsen 
	compacted = []
	for i in range(len(surfzunfold)):
		compacted.append([int(round(j/domainw)) for j in surfzunfold[i][0:2]]+[surfzunfold[i][2]])
	compacted = array(compacted)
	# group blocks
	mc = int(max(compacted[:,0])+1)
	nc = int(max(compacted[:,1])+1)
	averaged = [[[] for i in range(nc)] for j in range(mc)]
	for i in compacted:
		averaged[int(i[0])][int(i[1])].append(i[2])
	# average them
	for i in range(mc):
		for j in range(nc):
			averaged[i][j] = mean(averaged[i][j])
	# filter heights
	for i in range(mc):
		for j in range(nc):
			if (averaged[i][j] < 0 and heightfilter == 1) or (averaged[i][j] > 0 and heightfilter == -1):
				averaged[i][j] = 1
			else:
				averaged[i][j] = 0
	# scale proteins
	protcopy = list(array(list(prot)))
	for i in range(shape(protcopy)[0]):
		for j in range(2):
			protcopy[i][j] = int(round(protcopy[i][j]/domainw))
	protpts = [[0 for i in range(nc)] for j in range(mc)]
	protptsplus = [[0 for i in range(nc)] for j in range(mc)]
	# redundant frame border
	for i in protcopy:
		protpts[int(i[0])][int(i[1])] = 1
		for xd in range(-manhatdist,manhatdist+1):
			for yd in range(-manhatdist,manhatdist+1):
				protptsplus[int(i[0])+xd][int(i[1])+yd] = 1
	if 1:	
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
	write_framefilter(surfzunfold,array(domain),array(protptsplus),domainw,frameno)
	hinton_custom(array(averaged),array(protpts),array(protptsplus),array(domain),
		filename=(systemprefix+'-framefilter-figs/fig%05d.png'%frameno))

def batch_framefilter():
	for i in range(len(surfaces)):
		print 'Writing framefilter for frame: '+str(i)
		coarsen3(surfaces[i],proteins[i],6,i)

def view_average():
	mset.calculate_average_surface()
	meshplot(mset.surf_mean,vecs=mset.vecs[0])
	protpts = array([i-[0,0,25] for i in mean(mset.protein,axis=0)])
	meshpoints(protpts-[0,0,mean(mset.surf_position)],scale_factor=10,color=(1,1,1))

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#view_average()
batch_framefilter()


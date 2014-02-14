#!/usr/bin/env python

#---PLOTTERS: LIBRARIES, DEFINITIONS
#-------------------------------------------------------------------------------------------------------------

#---Plotting libraries
import numpy as np
import matplotlib as mpl

#---Remote operation if running batches
if plot_suppress:
	mpl.use('Agg')

#---Plotting libraries
#changed the following line to suit dark, which is having compatibility issues because new packages
import matplotlib.pyplot as plt
from matplotlib.pyplot import step
from matplotlib.ticker import NullFormatter
from matplotlib import rc
import pylab
from datetime import date
from mpl_toolkits.mplot3d import Axes3D

#---Import libraries
import numpy.polynomial
from mayavi import mlab

#---Color definitions
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors

#---Neutral fonts
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text', usetex=True) ########## this kills the program
#if you use amsmath, it gives you serif fonts and I don't know why argh
mpl.rc('text.latex', preamble='\usepackage{amsmath}')
#mpl.rc('text.latex', preamble='\usepackage{sfmath}')

#---PLOTTERS
#-------------------------------------------------------------------------------------------------------------

#---Mesh plotting functions via mayavi

def meshplot(data,vecs=None,translate=[0.0,0.0,0.0],show='both',wirecolor=None,
	tri=None,surfcolor=None,lsize=None,maxmin=None,opacity=1.):
	'''Plots a mesh given points in 3-space.'''
	if type(data) == list:
		data = array(data)
	if lsize == None: lsize=2.
	if maxmin == None:
		maxmin = [None,None]
	#---Automatically creates the mesh
	if tri == None:
		if shape(data)[1] != 3:
			data = unzipgrid(data,vecs=vecs)
		X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
		pts = mlab.points3d(X, Y, Z, Z)
		mesh = mlab.pipeline.delaunay2d(pts)
		pts.remove()
		if show == 'both':
			surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=(0,0,0),
				vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			surf2 = mlab.pipeline.surface(mesh,representation='surface',color=surfcolor,vmin=maxmin[0],
				vmax=maxmin[1],opacity=opacity,colormap='blue-red')
		elif show == 'wire':
			if wirecolor == None:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=(0,0,0),
					vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			elif wirecolor == 0:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=None,
					vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			else:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',
					line_width=lsize,color=wirecolor,vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,
					colormap='blue-red')
		elif show == 'surf':
			surf2 = mlab.pipeline.surface(mesh,representation='surface',color=surfcolor,line_width=lsize,
				vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
		mlab.axes(y_axis_visibility=False,x_axis_visibility=False,z_axis_visibility=False)
		mlab.xlabel("x")
		mlab.ylabel("y")
		mlab.zlabel("z")
	#---Custom mesh
	else:
		if shape(data)[1] != 3:
			data = unzipgrid(data,vecs=vecs)
		X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
		if show == 'both':
			surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',line_width=lsize,color=(0,0,0))
			surf2 = mlab.triangular_mesh(X,Y,Z,tri,representation='surface')
		elif show == 'wire':
			if wirecolor == None:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',
					line_width=lsize,color=(0,0,0))
			elif wirecolor == 0:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',line_width=lsize,color=None)
			else:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',
					line_width=lsize,color=wirecolor)
		elif show == 'surf':
			surf2 = mlab.triangular_mesh(X,Y,Z,tri,representation='surface')
		mlab.axes(y_axis_visibility=False,x_axis_visibility=False,z_axis_visibility=False)
		mlab.xlabel("x")
		mlab.ylabel("y")
		mlab.zlabel("z")

def meshpoints(data,vecs=None,translate=[0.0,0.0,0.0],scale_mode=None,scale_factor=None,color=(0,0,0),
	opacity=1,resolution=8):
	'''Add points to a mesh plot.'''
	if len(data) == 1 or (type(data) == numpy.ndarray and len(data) == 3):
		data = [data]
	if type(data) == list:
		data = array(data)
	X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
	if scale_factor == None:
		scale_factor = 3.
		scale_mode = 'none'
		s = Z
	elif type(scale_factor) == int:
		scale_factor = [scale_factor for i in range(len(data))]
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	else:
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	if shape(data)[1] != 3:
		data = unzipgrid(data,vecs=vecs)
	pts = mlab.points3d(X,Y,Z,s,scale_mode=scale_mode,mode='sphere',
		scale_factor=0.5,color=color,opacity=opacity,resolution=resolution)
		
def drawbox3d(vecs):
	'''Draw a bounding box.'''
	transpos = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
	seqs = [[0,1,2,3,0,4,5,1],[4,7,6,5],[7,3],[6,2]]
	for seq in seqs:
		dat = array([[transpos[s][i]*vecs[i] for i in range(3)] for s in seq])
		mlab.plot3d(dat[:,0],dat[:,1],dat[:,2],tube_radius=2, colormap='Spectral')

def unzipgrid(surf,vecs=None,grid=None,rounder_vecs=[1.0,1.0],reverse=0):
	'''Turns a 2D array into a set of points in 3-space. Duplicated in MembraneSet class.'''
	if type(surf) != ndarray:
		surf = array(surf)
	grid = [shape(surf)[i] for i in range(2)]
	if reverse != 0: grid = grid[::-1];
	if vecs != None:
		rounder_vecs = [vecs[i]/(grid[i]-1) for i in range(2)]
	replotin = surf
	surfreplot = []
	for i in range(grid[0]):
			for j in range(grid[1]):
				surfreplot.append([i*rounder_vecs[0],j*rounder_vecs[1],replotin[i,j]])
	surfreplot = array(surfreplot)
	return surfreplot

def meshplot_from_file(filename):
	'''Plots a triangulated mesh from a VTK-style file.'''
	data = mlab.pipeline.open(filename)
	surf = mlab.pipeline.surface(data)
	
def checkmesh(data,tri=None,vecs=None,grid=None,tess=1,show='both',wirecolor=None,surfcolor=None,lsize=None):
	'''Check a triangulated surface. Includes tesselation options.'''
	if type(data) != ndarray:
		data = array(data)
	if tess == 2:
		movetess = [[0,0],[0,1]]
	elif tess == -1:
		movetess = [[1,0]]
	elif tess == -1.5:
		movetess = [[-1,0]]
	elif tess == -2.5:
		movetess = [[-2,0]]
	elif tess == 3:
		movetess = [[0,0],[1,0],[0,1]]
	elif tess == -3:
		movetess = [[1,0],[0,1],[1,1]]
	elif tess == 4:
		movetess = [[0,0],[1,0],[0,1],[1,1]]
	elif tess == 9:
		movetess = [[0,0],[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]
	elif tess == 8:
		movetess = [[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]]
	elif tess == 'plus':
		movetess = [[1,0],[0,1],[-1,0],[0,-1]]
	else:
		movetess = [[0,0]]
	for move in movetess:
		meshplot(data,vecs=vecs,tri=tri,show=show,translate=[move[0]*vecs[0],move[1]*vecs[1],0],
			wirecolor=wirecolor,surfcolor=None,lsize=lsize)

#---Undulation plots

def plotter_undulations_2d_zoom(mset,logcolor=False,maxval=None,minval=None,
	imagefile=None,qmagfilter=[10**-6,10**6],lenscale=None):
	'''Single plot of the 2D undulation spectrum pattern, zoomed to the center.'''
	lenscale = mset.lenscale if lenscale == None else lenscale
	#---Plot the 2D spectrum
	[Lx,Ly] = [mean(mset.vecs[0])/lenscale,mean(mset.vecs[0])/lenscale]
	[m,n] = mset.griddims
	flanking = [int(round(i/2.)) for i in mset.griddims]
	widths = [i/8 for i in mset.griddims]
	fig = plt.figure(figsize=(5,5))
	ax = plt.subplot(111)
	ax.set_xlabel('x (nm)',fontsize=18)
	ax.set_ylabel('y (nm)',fontsize=18)
	ax.set_title('Undulation Spectrum (2D)',fontsize=18)
	x = [(j-m/2)/Lx*2*pi for j in range(0,n)]
	y = [(j-n/2)/Ly*2*pi for j in range(0,m)]
	X,Y = np.meshgrid(x,y)
	Z = mset.uqrawmean[int(flanking[0]-widths[0]):int(flanking[0]+widths[0]+1),
		int(flanking[1]-widths[1]):int(flanking[1]+widths[1]+1)]
	maxplotval = maxval if maxval != None else max([max(i) for i in mset.uqrawmean/2.])
	minplotval = minval if minval != None else 10**1
	if max([max(i) for i in mset.uqrawmean/2.]) > maxplotval: 
		print "Warning: your maximum color is lower than your true maximum! May be washing out the pattern."
	plt.imshow(Z,interpolation='nearest', cmap=mpl.cm.get_cmap('RdGy',100),
		norm=(mpl.colors.LogNorm(vmin=minplotval,vmax=maxplotval) if logcolor == True else None),
		origin='lower',aspect='equal',extent=[mset.griddims[0]/4,3/4*mset.griddims[0],
		mset.griddims[1]/4,3/4*mset.griddims[1]])
	plt.show()
	
def plotter_undulations_summary(mset,logcolor=False,maxval=None,minval=None,
	imagefile=None,qmagfilter=[10**-6,10**6],lenscale=None,zoom=False):
	'''Plot the 1D and 2D undulation spectra.'''
	ax2errors = False
	lenscale = mset.lenscale if lenscale == None else lenscale
	#---Plot the 2D spectrum
	[Lx,Ly] = [mean(mset.vecs[0])/lenscale,mean(mset.vecs[0])/lenscale]
	[m,n] = mset.griddims
	fig = plt.figure(figsize=(15, 5))
	ax = plt.subplot2grid((1,2+(1 if zoom==True else 0)),(0,1))
	ax.set_xlabel('x (nm)',fontsize=18)
	ax.set_ylabel('y (nm)',fontsize=18)
	ax.set_title('Undulation Spectrum (2D)',fontsize=18)
	x = [(j-m/2)/Lx*2*pi for j in range(0,n)]
	y = [(j-n/2)/Ly*2*pi for j in range(0,m)]
	X,Y = np.meshgrid(x,y)
	Z = mean(mset.uqcollect,axis=0)
	maxplotval = maxval if maxval != None else max(mset.uqrawmean)
	minplotval = minval if minval != None else 10**1
	if max(mset.uqrawmean) > maxplotval: 
		print "Warning: your maximum color is lower than your true maximum! May be washing out the pattern."
	plt.imshow(Z,interpolation='nearest', cmap=mpl.cm.get_cmap('RdGy',100),
		norm=(mpl.colors.LogNorm(vmin=minplotval,vmax=maxplotval) if logcolor == True else None),
		origin='lower',aspect='equal',extent=[0,mset.griddims[0],0,mset.griddims[1]])
	xtickx = [int(i) for i in list(arange(0,mset.griddims[0],int(mset.griddims[0]/10)))]
	xticky = [int(i/lenscale) for i in list(arange(0,mset.vecs[0][0],int(mset.vecs[0][0]/(mset.griddims[0]-1)*
		int(mset.griddims[0]/10))))]
	ytickx = [int(i) for i in list(arange(0,mset.griddims[1],int(mset.griddims[1]/10)))]
	yticky = [int(i/lenscale) for i in list(arange(0,mset.vecs[0][1],int(mset.vecs[0][1]/(mset.griddims[1]-1)*
		int(mset.griddims[1]/10))))]
	plt.xticks(xtickx,xticky)
	plt.yticks(ytickx,yticky)
	#---Plot the 1D spectrum
	ax2 = plt.subplot2grid((1,2+(1 if zoom==True else 0)), (0,0))
	ax2.set_xlabel(r"$\left|q\right|$",fontsize=18)
	ax2.set_ylabel(r"$\left\langle z_{q}z_{-q}\right\rangle$",fontsize=18)
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.grid(True)
	#plt.title('Undulation Spectrum (1D)',fontsize=18)
	spectrum = array(sorted(zip(*[mset.qrawmean,mset.uqrawmean,mset.uqrawstd]), key=lambda x: x[0]))[1:]
	if ax2errors:
		ax2.errorbar(spectrum[:,0],spectrum[:,1],yerr=spectrum[:,2],color='b',fmt='.',alpha=0.2)
	ax2.scatter(spectrum[:,0],spectrum[:,1],marker='o',color='k')
	#---Fitting
	spectrumf = array(filter(lambda x: x[0] >= qmagfilter[0] and x[0] <= qmagfilter[1], spectrum))
	spectrumf2 = array(filter(lambda x: x[0] >= 1./10.*qmagfilter[0] and x[0] <= qmagfilter[1]*10., spectrum))
	[bz,az]=numpy.polyfit(log(spectrumf[:,0]),log(spectrumf[:,1]),1)
	#---Compute and display the bending rigidity
	area = double(mean([i[0]*i[1] for i in mset.vecs])/lenscale**2)
	print 'q-magnitude scaling = '+str(bz)
	kappa = 1/exp(az)/area
	print 'kappa = '+str(kappa)
	leftcom = [mean(log(spectrumf[:,0])),mean(log(spectrumf[:,1]))]
	az_enforced = leftcom[1]+4.*leftcom[0]
	kappa_enforced = 1./exp(az_enforced)/area
	print 'kappa_enforced = '+str(kappa_enforced)
	#---Fitting
	ymod=[exp(az)*(i**bz) for i in spectrumf[:,0]]
	xmod=[i for i in spectrumf[:,0]]
	ymod2=[exp(az)*(i**bz) for i in spectrumf2[:,0]]
	xmod2=[i for i in spectrumf2[:,0]]
	ymod3=[exp(az)*(i**bz) for i in spectrumf2[:,0]]
	xmod3=[i for i in spectrumf2[:,0]]
	ymod4=[exp(az_enforced)*(i**-4.) for i in spectrumf2[:,0]]
	xmod4=[i for i in spectrumf2[:,0]]
	#---Plot fitted lines
	ax2.plot(xmod2,ymod2,color='#FF3399',linestyle='-',linewidth=2.5)
	ax2.plot(xmod4,ymod4,color='#3399FF',linestyle='-',linewidth=2.5)
	#ax2.plot(xmod,ymod,marker='.',color='w')
	#---Write bending rigidity on the plot
	#ax2.text(0.1, 0.2, r'$\kappa_{'+str('%1.2f'%bz)+'} = %3.2f$'%kappa, transform=ax2.transAxes,fontsize=16)
	ax2.text(0.1, 0.2, r'$\kappa = %3.2f$'%kappa, transform=ax2.transAxes,fontsize=16)
	#ax2.text(0.1, 0.1, r"$\kappa_{-4} = %3.2f$"%kappa_enforced, transform=ax2.transAxes,fontsize=16)
	#---Save and show
	plt.tight_layout()
	if imagefile != None:
		plt.savefig(imagefile,dpi=500)
	if zoom == True:
		lenscale = mset.lenscale if lenscale == None else lenscale
		#---Plot the 2D spectrum, zoomed in to the center pattern
		[Lx,Ly] = [mean(mset.vecs[0])/lenscale,mean(mset.vecs[0])/lenscale]
		[m,n] = mset.griddims
		flanking = [int(round(i/2.))-1 for i in mset.griddims]
		widths = [i/6 for i in mset.griddims]
		ax3 = plt.subplot2grid((1,3),(0,2))
		ax3.set_xlabel('x',fontsize=18)
		ax3.set_ylabel('y',fontsize=18)
		ax3.set_title('Undulation Spectrum (2D)',fontsize=18)
		x = [(j-m/2)/Lx*2*pi for j in range(0,n)]
		y = [(j-n/2)/Ly*2*pi for j in range(0,m)]
		X,Y = np.meshgrid(x,y)
		Z = mean(mset.uqcollect,axis=0)[int(flanking[0]-widths[0]):int(flanking[0]+widths[0])+1,
			int(flanking[1]-widths[1]):int(flanking[1]+widths[1])+1]
		maxplotval = maxval if maxval != None else max(mset.uqrawmean)
		minplotval = minval if minval != None else 10**1
		if max(mset.uqrawmean) > maxplotval:
			print "Warning: your maximum color is lower than your true maximum! Pattern may be washed-out."
		plt.imshow(Z,interpolation='nearest', cmap=mpl.cm.get_cmap('RdGy',100),
			norm=(mpl.colors.LogNorm(vmin=minplotval,vmax=maxplotval) if logcolor == True else None),
			origin='lower',aspect='equal',extent=[0,shape(Z)[0],0,shape(Z)[1]])
		plt.xticks(range(shape(Z)[0]),[(i+flanking[0]-1) for i in range(shape(Z)[0])])
		plt.yticks(range(shape(Z)[1]),[(i+flanking[1]-1) for i in range(shape(Z)[1])])
	plt.show()

#---Hinton plots
	
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
				color = '#FFFF33' ##### yellow
			elif w == 1:
				color = '#0066CC' #blue 
			hinton_blob(_x - 0.5, height - _y + 0.5, 1,color)
			if domain[y,x] == 1:
				hinton_blob(_x - 0.5, height - _y + 0.5, 1,'#FF9900') # orange
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
		
		
		
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------		
		

	
#---PLOTTING FUNCTIONS, HACKED
#-------------------------------------------------------------------------------------------------------------
	
#---Radial distribution function plots

############# MIRROR VERSION - KEEP THIS IN THE CODE AND CLEAN IT UP
#def plot_gr_voronoi(self,label,specs=None,c=0,savefig=None,style=None,filename=None,direction='z'):
if 0:
	self=mset
	specs=None
	c=0
	savefig=None
	style=None
	filename=None
	direction='z'
	label0='DOPC P-to-MG'
	label1='DOPC P-to-CL'
	binwidth = 2.
	cutoff = mean(mset.vecs,axis=0)[2]/2.
	spec0 = ['direction','z']
	specs = spec0+specs if specs != None else spec0
	pdist0 = self.getdata(label0).get(specs)
	pdist1 = self.getdata(label1).get(specs)
	hist0 = numpy.histogram(pdist0,range=(0,cutoff),bins=int(cutoff/binwidth+1))
	hist1 = numpy.histogram(pdist1,range=(0,cutoff),bins=int(cutoff/binwidth+1))
	plt.bar(hist0[1][:-1],2*hist0[0],color='b',alpha=0.5)
	plt.bar(hist1[1][:-1],-1*hist1[0],color='r',alpha=0.5)
	plt.plot(hist1[1][:-1],2*hist0[0]-hist1[0],'o-')
	plt.xlim(0,cutoff)
	plt.show()

########### CLEAN THIS UP AND USE AS STANDARD gr-voronoi HISTOGRAM FUNCTION
#def plot_gr_voronoi(self,label,specs=None,c=0,savefig=None,style=None,filename=None,direction='z'):
if 0:
	self=mset
	label=None
	specs=None
	c=0
	savefig=None
	style=None
	filename=None
	direction='z'
	'''Plots a histogram of data from the gr-voronoi calculation.'''
	if label == None: label = self.avail()
	if type(label) == list and type(label[0]) == int:
		label = [self.avail()[i] for i in label]
	if type(label) == str: label = [label]
	spec0 = ['direction','z']
	specs = spec0+specs if specs != None else spec0
	print specs
	#specs = ['direction',dircode].extend(specs) if specs != None else ['direction',dircode]
	fig, ax = plt.subplots()
	histplots = []
	print specs
	for i in range(len(label)):
		name = label[i]
		c = clrs[i%len(clrs)]
		c = clrs[1]
		pdistsz = self.getdata(name).get(specs)
		print pdistsz
		cutoff = mean(mset.vecs,axis=0)[2]/2.
		if style == None:
			thisplot = hist(pdistsz, bins=40,range=(0,cutoff),color=c, 
				normed=True, alpha=(0.5 if len(label) > 1 else None),label=name)
		else:
			thisplot = hist(pdistsz, bins=40,range=(0,cutoff),color=c, 
				normed=True,histtype='step',linewidth=3,label=name)
		histplots.append(thisplot)
	ax.set_xlim((0,cutoff))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Distance (\AA)', labelpad = 10)
	plt.ylabel('Number Frequency', labelpad = 10)
	plt.title('RDF of Lipid-Ion distances, Voronoi Bins')
	fig.tight_layout()
	plt.legend()
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=300)
		else:
			plt.savefig(filename+'.png', dpi=300)
	plt.show()

	
#---Radial distribution function plots

def plot_gr_voronoi2(self,label,specs=None,c=0,savefig=None,style=None,filename=None,direction='z'):
	'''Plots a histogram of data from the gr-voronoi calculation.'''
	if label == None: label = self.avail()
	if type(label) == list and type(label[0]) == int:
		label = [self.avail()[i] for i in label]
	if type(label) == str: label = [label]
	if direction == 'z':
		dircode = 0
	elif direction == 'lateral':
		dircode = 1
	specs = ['direction',dircode].extend(specs) if specs != None else ['direction',dircode]
	fig, ax = plt.subplots()
	histplots = []
	for i in range(len(label)):
		name = label[i]
		c = clrs[i%len(clrs)]
		print specs
		pdistsz = self.getdata(name).get([specs]) ##### bracket hack
		if style == None:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=True, alpha=(0.5 if len(label) > 1 else None),label=name)
		else:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=True,histtype='step',linewidth=3,label=name)
		histplots.append(thisplot)
	ax.set_xlim((0, 60))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Distance (\AA)', labelpad = 10)
	plt.ylabel('Number Frequency', labelpad = 10)
	plt.title('RDF of Lipid-Ion distances, Voronoi Bins')
	fig.tight_layout()
	plt.legend()
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=300)
		else:
			plt.savefig(filename+'.png', dpi=300)
	plt.show()

	
def plot_gr_voronoi_heat(self,label,specs=None,filename=None,savefig=False,logcolor=True):
	'''Heat map of lateral vs normal distances from gr-voronoi calculation.'''
	if type(label) != str:
		label = mset.avail()[label]
	fig = plt.figure()
	farcutoff = 1./2*mean(array(self.vecs)[:,2])-20.
	basespec = ['direction',0]
	specsz = [basespec] + [specs] if specs != None else basespec
	pdistsz = self.getdata(label).get(specsz)
	pdistsz_mask = np.ma.array (pdistsz, mask=np.isnan(pdistsz))
	basespec = ['direction',1]
	specs2d = [basespec] + [specs] if specs != None else basespec        
	########### check that you don't care how many lists on the line above...
	print specs2d
	pdists2d = self.getdata(label).get(specs2d)
	pdists2d_mask = np.ma.array (pdists2d, mask=np.isnan(pdists2d))
	H, xedges, yedges = histogram2d(pdistsz_mask,pdists2d_mask,bins=50)
	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	minplotval = 1.0
	maxplotval = max([max(pdistsz),max(pdists2d)])
	cmap = mpl.cm.jet
	cmap.set_bad('k',1.)
	plt.imshow(H, extent=None, interpolation='nearest',aspect='equal',origin='lower',
		norm=(mpl.colors.LogNorm(vmin=minplotval,vmax=maxplotval) if logcolor == True else None),cmap=cmap)
	plt.colorbar()
	ax = fig.add_subplot(111)
	ax.set_xlabel('Lateral distance (\AA)',fontsize=16)
	ax.set_ylabel('Normal distace (\AA)',fontsize=16)
	ax.set_title('RDF Histogram, Voronoi Bins')
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=100)
		else:
			plt.savefig(filename+'.png', dpi=100)
	plt.show()
	
def plot_aplxxxx(self,label='triangles',ltypes=None,specs=None,filename=None,savefig=False,monos=[0,1]):
	'''Area per lipid calculations and plot.'''
	if type(label) != str:
		label = self.avail()[label]
	lipid_areas = []
	if ltypes == None:
		ltypes = range(len(self.resnames))
	for ltype in ltypes:
		areas = []
		for mono in monos:
			print 'Calculating areas: lipid = '+self.resnames[ltype]+' monolayer '+str(mono)+'.'
			for fr in range(len(self.getdata('triangles').label)):
				points = array(self.getdata(label).get([['type','points'],['frame',fr],['monolayer',mono]]))
				tri1 = array(self.getdata(label).get([['type','lines'],['frame',fr],['monolayer',0]]))
				validres = [self.monolayer_residues[mono].index(i) 
					for i in self.monolayer_by_resid[mono][ltype]]
				areas.append([sum([1./6*norm(cross(t[0][0:2]-t[1][0:2],t[2][0:2]-t[1][0:2])) 
					for t in [[points[j] for j in tri1[i]] for i in list(where(any(tri1==p,axis=1))[0])]]) 
					for p in validres])
		lipid_areas.append(areas)
	areas = array(lipid_areas)
	fig, ax = plt.subplots()
	histplots = []
	for d in range(len(lipid_areas)):
		c = clrs[d%len(clrs)]
		thisplot = hist(mean(lipid_areas[d],axis=0),normed=True,label=self.resnames[d],alpha=0.5,color=c)
		histplots.append(thisplot)
		ax.text(0.7,0.6-d*0.05, r"$A_{"+str(self.resnames[d])+"}= "+str('%2.2f'%mean(lipid_areas[d]))+
			"\pm"+str('%1.2f'%std(lipid_areas[d]))+"$",transform=ax.transAxes,fontsize=12)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Area per lipid ($\AA^{2}$)',labelpad = 10,fontsize=16)
	plt.ylabel('Normed frequency',labelpad = 10,fontsize=16)
	plt.title('Area per lipid',fontsize=16)
	fig.tight_layout()
	plt.legend()
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=100)
		else:
			plt.savefig(filename+'.png', dpi=100)
	plt.show()
	
def plot_apl_compare(self0,self1,ltypes1=None,ltypes0=None,label='triangles',specs=None,
	filename=None,savefig=False,monos=[0,1],name0='',name1='',gaussfit=True):
	'''Area per lipid calculations and plot to compare two simulations.'''
	#---Simulation 0 area calculation
	if ltypes0 == None:
		ltypes0 = range(len(self0.resnames))
	elif type(ltypes0) == str:
		ltypes0 = [self0.resnames.index(ltypes0)]
	if type(label) != str:
		label = self0.avail()[label]
	lipid_areas0 = []
	for ltype in ltypes0:
		areas0 = []
		for mono in monos:
			print 'Calculating areas: lipid = '+self0.resnames[ltype]+' monolayer '+str(mono)+'.'
			for fr in range(len(self0.getdata('triangles').label)):
				points = array(self0.getdata(label).get([['type','points'],['frame',fr],['monolayer',mono]]))
				tri1 = array(self0.getdata(label).get([['type','lines'],['frame',fr],['monolayer',mono]]))
				validres = [self0.monolayer_residues[mono].index(i) 
					for i in self0.monolayer_by_resid[mono][ltype]]
				areas0.append([sum([1./6*norm(cross(t[0][0:2]-t[1][0:2],t[2][0:2]-t[1][0:2])) 
					for t in [[points[j] for j in tri1[i]] for i in list(where(any(tri1==p,axis=1))[0])]]) 
					for p in validres])
		lipid_areas0.append(areas0)
	areas0 = array(lipid_areas0)
	#---Simulation 1 area calculation
	if ltypes1 == None:
		ltypes1 = range(len(self1.resnames))
	elif type(ltypes1) == str:
		ltypes1 = [self1.resnames.index(ltypes1)]
	if type(label) != str:
		label = self1.avail()[label]
	lipid_areas1 = []
	for ltype in ltypes1:
		areas1 = []
		for mono in monos:
			print 'Calculating areas: lipid = '+self1.resnames[ltype]+' monolayer '+str(mono)+'.'
			for fr in range(len(self1.getdata('triangles').label)):
				points = array(self1.getdata(label).get([['type','points'],['frame',fr],['monolayer',mono]]))
				tri1 = array(self1.getdata(label).get([['type','lines'],['frame',fr],['monolayer',mono]]))
				validres = [self1.monolayer_residues[mono].index(i) 
					for i in self1.monolayer_by_resid[mono][ltype]]
				areas1.append([sum([1./6*norm(cross(t[0][0:2]-t[1][0:2],t[2][0:2]-t[1][0:2])) 
					for t in [[points[j] for j in tri1[i]] for i in list(where(any(tri1==p,axis=1))[0])]]) 
					for p in validres])
		lipid_areas1.append(areas1)
	areas1 = array(lipid_areas1)
	#---Combination plots
	fig, ax = plt.subplots()
	histplots = []
	for d in range(len(lipid_areas0)):
		c = clrs[d%len(clrs)]
		thisplot = hist(mean(lipid_areas0[d],axis=0),normed=True,label=name0+' '+self0.resnames[ltypes0[d]],
			alpha=0.5,color=c)
		histplots.append(thisplot)
		ax.text(0.7,0.6-d*0.05, r"$A_{"+str(self0.resnames[ltypes0[d]])+"}= "+
			str('%2.2f'%mean(lipid_areas0[d]))+"\pm"+str('%1.2f'%std(lipid_areas0[d]))+"$",
			transform=ax.transAxes,fontsize=12)
		y = normpdf(thisplot[1], mean(lipid_areas0[d],axis=0), std(lipid_areas0[d],axis=0))
		plt.plot(bins, y, 'r--')

	for d in range(len(lipid_areas1)):
		c = clrs[(d+len(ltypes0))%len(clrs)]
		thisplot = hist(mean(lipid_areas1[d],axis=0),normed=True,label=name1+' '+self1.resnames[ltypes1[d]],
			alpha=0.5,color=c)
		histplots.append(thisplot)
		ax.text(0.7,0.6-(d+len(ltypes0))*0.05, r"$A_{"+str(self1.resnames[ltypes1[d]])+"}= "+
			str('%2.2f'%mean(lipid_areas1[d]))+"\pm"+str('%1.2f'%std(lipid_areas1[d]))+"$",
			transform=ax.transAxes,fontsize=12)
		y = normpdf(thisplot[1], mean(lipid_areas0[d],axis=0), std(lipid_areas0[d],axis=0))
		plt.plot(bins, y, 'r--')

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Area per lipid ($\AA^{2}$)',labelpad = 10,fontsize=16)
	plt.ylabel('Normed frequency',labelpad = 10,fontsize=16)
	plt.title('Area per lipid',fontsize=16)
	fig.tight_layout()
	plt.legend()
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=100)
		else:
			plt.savefig(filename+'.png', dpi=100)
	plt.show()
	
def plot_apl_compare_stacked(selfs,ltypes,label='triangles',specs=None,
	filename=None,savefig=False,monos=[0,1],names=None,gaussfit=True):
	'''Area per lipid calculations and plot to compare two simulations.'''
	all_areas = []
	for self0 in selfs:
		#---Simulation 0 area calculation
		lipid_areas0 = []
		for ltype in ltypes:
			areas0 = []
			for mono in monos:
				print 'Calculating areas: lipid = '+self0.resnames[ltype]+' monolayer '+str(mono)+'.'
				for fr in range(len(self0.getdata('triangles').label)):
					points = array(self0.getdata(label).get([['type','points'],['frame',fr],
						['monolayer',mono]]))
					tri1 = array(self0.getdata(label).get([['type','lines'],['frame',fr],
						['monolayer',mono]]))
					validres = [self0.monolayer_residues[mono].index(i) 
						for i in self0.monolayer_by_resid[mono][ltype]]
					areas0.append([sum([1./6*norm(cross(t[0][0:2]-t[1][0:2],t[2][0:2]-t[1][0:2])) 
						for t in [[points[j] for j in tri1[i]] for i in list(where(any(tri1==p,axis=1))[0])]]) 
						for p in validres])
			lipid_areas0.append(areas0)
		allareas.append(lipid_areas0)
		
	if 0:
		#---Combination plots
		fig, ax = plt.subplots()
		histplots = []
		for d in range(len(lipid_areas0)):
			c = clrs[d%len(clrs)]
			thisplot = hist(mean(lipid_areas0[d],axis=0),normed=True,label=name0+' '+
				self0.resnames[ltypes0[d]],alpha=0.5,color=c)
			histplots.append(thisplot)
			ax.text(0.7,0.6-d*0.05, r"$A_{"+str(self0.resnames[ltypes0[d]])+"}= "+
				str('%2.2f'%mean(lipid_areas0[d]))+"\pm"+str('%1.2f'%std(lipid_areas0[d]))+"$",
				transform=ax.transAxes,fontsize=12)
			y = normpdf(thisplot[1], mean(lipid_areas0[d],axis=0), std(lipid_areas0[d],axis=0))
			plt.plot(bins, y, 'r--')

		for d in range(len(lipid_areas1)):
			c = clrs[(d+len(ltypes0))%len(clrs)]
			thisplot = hist(mean(lipid_areas1[d],axis=0),normed=True,label=name1+' '+
				self1.resnames[ltypes1[d]],alpha=0.5,color=c)
			histplots.append(thisplot)
			ax.text(0.7,0.6-(d+len(ltypes0))*0.05, r"$A_{"+str(self1.resnames[ltypes1[d]])+"}= "+
				str('%2.2f'%mean(lipid_areas1[d]))+"\pm"+str('%1.2f'%std(lipid_areas1[d]))+"$",
				transform=ax.transAxes,fontsize=12)
			y = normpdf(thisplot[1], mean(lipid_areas0[d],axis=0), std(lipid_areas0[d],axis=0))
			plt.plot(bins, y, 'r--')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_position(('outward', 20))
		ax.spines['left'].set_position(('outward', 30))
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		plt.xlabel('Area per lipid ($\AA^{2}$)',labelpad = 10,fontsize=16)
		plt.ylabel('Normed frequency',labelpad = 10,fontsize=16)
		plt.title('Area per lipid',fontsize=16)
		fig.tight_layout()
		plt.legend()
		if savefig: 
			if not filename:
				plt.savefig(name+'.png', dpi=100)
			else:
				plt.savefig(filename+'.png', dpi=100)
		plt.show()
	
def plot_gr_voronoi_compare(self0,self1,label0=None,label1=None,specs=None,colors=None,savefig=None,
	style=None,filename=None,direction='z'):
	'''Plots a histogram of data from the gr-voronoi calculation.'''

	#---Simulation 0
	if direction == 'z':
		dircode = 0
	elif direction == 'lateral':
		dircode = 1
	specs = ['direction',dircode].extend(specs) if specs != None else ['direction',dircode]
	fig, ax = plt.subplots()
	histplots = []
	for i in range(len(label0)):
		name = label0[i]
		c = clrs[colors[0]]
		print specs
		pdistsz = self0.getdata(name).get(specs)
		if style == None:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=False,weights=[1./len(self0.vecsi) for i in range(len(pdistsz))],
				alpha=(0.5),label=name)
			print [sum(thisplot[0][:i]) for i in range(len(thisplot[0])-1)]
		else:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=False,weights=[1./len(self0.vecsi) for i in range(len(pdistsz))],
				histtype='step',linewidth=3,label=name)
			print [sum(thisplot[0][:i]) for i in range(len(thisplot[0])-1)]
		histplots.append(thisplot)

	#---Simulation 1
	if direction == 'z':
		dircode = 0
	elif direction == 'lateral':
		dircode = 1
	specs = ['direction','z']
	print label1
	for i in range(len(label1)):
		name = label1[i]
		c = clrs[colors[1]]
		print specs
		pdistsz = self1.getdata(name).get(specs)
		if style == None:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=False, weights=[1./len(self0.vecsi) for i in range(len(pdistsz))],alpha=(0.5),
				label=name)
			print [sum(thisplot[0][:i]) for i in range(len(thisplot[0])-1)]
		else:
			thisplot = hist(pdistsz, bins=40, range=(0,60), color=c, 
				normed=False,weights=[1./len(self0.vecsi) for i in range(len(pdistsz))],histtype='step',
				linewidth=3,label=name)
			print [sum(thisplot[0][:i]) for i in range(len(thisplot[0])-1)]
		histplots.append(thisplot)


	ax.set_xlim((0, 60))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Distance (\AA)', labelpad = 10)
	#plt.ylabel('Number Frequency', labelpad = 10)
	plt.ylabel('Counts', labelpad = 10)
	plt.title('RDF of Lipid-Ion distances, Voronoi Bins')
	fig.tight_layout()
	plt.legend()
	
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=300)
		else:
			plt.savefig(filename+'.png', dpi=300)
	#plt.show()
	
	return histplots

	
#---UNFINISHED PLOTTERS BEWARE
#-------------------------------------------------------------------------------------------------------------

#---Radial distribution function plots
	
def plotter_gr(mset,which_data=None,outfile=None,grresolution=2.0):
	'''UNFINISHED.'''
	if which_data == None:
		which_data = range(len(mset.gr_data))
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel('g(r) (A)')
	ax.set_ylabel('Density')
	ax.set_title('Lipid-Ion g(r)')
	for w in range(len(which_data)):
		farcutoff = 1./2*mean(mset.vecs[:,2])-20.
		pair_distances = concatenate([i.reshape(-1) for i in mset.gr_data[which_data[w]]])
		hist, bin_edges = histogram(pair_distances, bins=arange(0,farcutoff,grresolution))
		mid = (bin_edges[1:] + bin_edges[:-1])/2
		#histnorm = [hist[i]/(4.*pi*((mid[i])**2)) for i in range(len(hist))]
		c = cm.Dark2(w/float(len(which_data)),1)
		plt.plot(mid, hist, 'o-',linewidth=2,label=mset.gr_data_labels[which_data[w]],color=c)
	plt.legend(tuple(mset.gr_data_labels))
	if outfile != None:
		savefig(outfile)
	plt.show()
	
def plot_gr_v2(self,grresolution=2.0):
	'''UNFINISHED.'''
	defcolors = ['blue','purple','black']
	pairdists = self.gr_pair_distances
	boxvecs = self.universe.trajectory[0].dimensions
	boxv = array([boxvecs])
	meandens = len(pairdists[0])/mean(boxv[:,0])/mean(boxv[:,1])
	farcutoff = 1/2.*(mean(boxv[:,0])+mean(boxv[:,1]))/2.
	subject_by_frame = []
	for framedists in pairdists:
		subject_by_frame.append(array(filter(lambda x: x != 0.0,[i for j in framedists for i in j])))
	subject = array([k for l in subject_by_frame for k in l])
	hist, bin_edges = histogram(subject, bins=arange(0,farcutoff,grresolution))
	areas = [pi*((k+resolution)**2-k**2) for k in bin_edges[:-1]]
	histnorm = [normfac*hist[i]/0.5/areas[i]/meandens/len(pairdists[0])/len(pairdists) \
		for i in range(len(hist))]
	mid = (bin_edges[1:] + bin_edges[:-1])/2
	plt.plot(mid, histnorm, 'o-',color='black',linewidth=2,)
	plt.show()
	
def plot_gr_voronoi_deprecated(self,which_data=None):
	'''UNFINISHED.'''
	if which_data == None:
		which_data = -1
	pdistsz = array(self.gr_data[which_data][0]).flatten()
	pdists2d = array(self.gr_data[which_data][1]).flatten()
	H, xedges, yedges = histogram2d([i for j in pdistsz for i in j], [i for j in pdists2d for i in j], \
		bins=40)
	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	plt.imshow(H, extent=None, interpolation='nearest',aspect='equal')
	plt.colorbar()
	plt.show()

def plotter_gr_dual(mset,which_data=None,outfile=None,grresolution=2.0):
	'''UNFINISHED.'''
	if which_data == None:
		which_data = range(len(mset.gr_data))
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel('g(r) group1 (A)')
	ax.set_ylabel('g(r) group2 (A)')
	ax.set_title('Lipid-Ion g(r)')
	for w in which_data[0:1]:
		farcutoff = 1./2*mean(mset.vecs[:,2])-20.
		pair_distances = concatenate([i.reshape(-1) for i in mset.gr_data[w]])
		pair_distances2 = concatenate([i.reshape(-1) for i in mset.gr_data[w+1]])
		hist, xedges, yedges = histogram2d(pair_distances, pair_distances2, \
			bins=arange(0,farcutoff+grresolution,grresolution))
		histnorm = [[hist[i][j]/(4.*pi*(xedges[i+1]**2+yedges[i+1]**2)) for i in range(len(hist))] \
			for j in range(len(hist[0]))]
		extent = [xedges[1], xedges[-1], yedges[-1], yedges[1]]
		imshow(histnorm, extent=extent, interpolation='nearest',origin='lower')
		plt.colorbar()
		if outfile != None:
			savefig(outfile)
		plt.show()
		
#---Basic plots

def plotter_histogram(inputdat):
	'''UNDER CONSTRUCTION.'''
	fig = plt.figure()
	ax2 = fig.add_subplot(111)
	#for l in range(0,len(inputdat)):
	#	subject = array(inputdat[l])  	
	#	n, bins, patches = ax2.hist(subject,range=(0,max(subject)),normed=1,alpha=0.4)
	subject = array(inputdat)  	
	n, bins, patches = ax2.hist(subject,range=(min(subject),max(subject)),alpha=0.4,bins=10)
	ax2.set_title('title')
	ax2.set_xlabel('xlabel')
	ax2.set_ylabel('ylabel')
	ax2.grid(True)
	plt.show()
	# use this for a standalone histogram: n, bins, patches = hist(tmp, 50, normed=1, histtype='stepfilled') 

#---Simple 3D plotting functions from matplotlib

def plot3d(surf,movezeros=0,zrange=0,plot2d=0):
	'''UNDER CONSTRUCTION.'''
	print 'in plot3d'
	if type(surf) is list:	
		surf = array(surf)
	m = shape(surf)[0]
	n = shape(surf)[1]
	if n == 3:
		plot3dstuff = surf
	else:
		replotin = surf
		surfreplot = []
		for i in range(m):
			    for j in range(n):
					surfreplot.append([i,j,replotin[i,j]])
		surfreplot = array(surfreplot)
		if movezeros != 0:
			plot3dstuff = array([i for i in surfreplot if i[2] != 0.0])
		else:
			plot3dstuff = surfreplot
	if plot2d == 0:
		fig = plt.figure()
		ax = fig.add_subplot(111,projection='3d')
		#ax.set_xlim3d([0,infty]);
		#if zrange != 0:
		#	ax.set_zlim3d(zrange[0],zrange[1])
		ax.scatter(plot3dstuff[:,0],plot3dstuff[:,1],plot3dstuff[:,2],marker='x')
		plt.show()
	else:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.scatter(plot3dstuff[:,0],plot3dstuff[:,1],marker='x')
		plt.show()

def plotter_contour(inputdat):
	'''UNDER CONSTRUCTION.'''
	m = shape(inputdat)[0]
	n = shape(inputdat)[1]
	replotin = array(inputdat)
	Z=[]
	for i in range(n):
		tempz = []
		for j in range(m):
			tempz.append(replotin[j,i])
		Z.append(tempz)
	Z=array(Z)
	X=array([[i for i in range(m)] for j in range(n)])	#!!!!!!!!!!!!!!!
	Y=array([[j for i in range(m)] for j in range(n)])
	fig = plt.figure(figsize=(8*1.2, 6*1.2),)
	ax = fig.gca()
	cset = ax.contour(X, Y, Z,cmap=cm.Set2)
	im = imshow(Z,cmap=cm.coolwarm)
	plt.colorbar(im, orientation='horizontal', shrink=0.5)
	ax.set_xlabel('X (nm)')
	ax.set_xlim(0,m-1)
	ax.set_ylabel('Y (nm)')
	ax.set_ylim(0,n-1)
	plt.show()

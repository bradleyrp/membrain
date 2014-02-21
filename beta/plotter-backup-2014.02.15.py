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
#---Nb: changed the following line to suit dark, which is having compatibility issues because new packages
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
		
def unzipgrid(surf,vecs=None,grid=None,rounder_vecs=[1.0,1.0],reverse=0):
	'''Turns a 2D array into a set of points in 3-space for meshplot and meshpoints.'''
	#---Nb: this function is also a child of the MembraneSet class.
	#---Nb: duplicated here for plotting purposes only.
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
		
def drawbox3d(vecs):
	'''Draw a bounding box.'''
	transpos = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
	seqs = [[0,1,2,3,0,4,5,1],[4,7,6,5],[7,3],[6,2]]
	for seq in seqs:
		dat = array([[transpos[s][i]*vecs[i] for i in range(3)] for s in seq])
		mlab.plot3d(dat[:,0],dat[:,1],dat[:,2],tube_radius=2, colormap='Spectral')

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
	print 'lenscale = '+str(lenscale)
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


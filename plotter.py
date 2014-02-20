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

#---Use the axesgrid toolkit for insets
import mpl_toolkits.axes_grid.inset_locator

#---Tools for importing VTU files
import numpy, glob, sys, scipy, os
import xml.etree.ElementTree as ET
from numpy.linalg import norm

#---Gridspec tools
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

#---undulation plots

def plotter2d(ax,mset,dat=None,nlabels=None,tickshow=False,cmap=None,lims=None,
	cmap_washout=None,lognorm=True,inset=True,i2wid=1,ticklabel_show=[1,1],label_style=None):
	'''Standard function for plotting 2D spectra from MembraneSet objects.'''
	if dat == None: dat = mset.undulate_hqhq2d
	if cmap == None: cmap=mpl.cm.binary; cmap_washout = 0.65;
	if cmap_washout == None: cmap_washout = 1.0
	if lims == None: lims = [None,None]
	if nlabels == None: nabels = 3
	if tickshow == True: tickshow = [1,1]
	if label_style == 'xy':
		xlabel = r'${x(\mathrm{nm})$'
		ylabel = r'${y(\mathrm{nm})$'
	elif label_style == 'q':
		xlabel = r'${\pi\left|q_x\right|}^{-1}(\mathrm{nm})$'
		ylabel = r'${\pi\left|q_y\right|}^{-1}(\mathrm{nm})$'
	else:
		ticklabel_show = [0,0]		
	lognorm = mpl.colors.LogNorm() if lognorm == True else None
	m,n = shape(dat)
	#---follows scipy DFFT convention on even/odd location of Nyquist component
	cm,cn = [int(i/2) for i in shape(dat)]
	#---nlabels sets the number of printed labels
	if (m < 5 or n < 5): nlabels = 1; sskipx = 1;sskipy =1
	elif nlabels == None: 
		nlabels = 3
		sskipx = m/nlabels
		sskipy = n/nlabels
	else:
		sskipx = m/nlabels
		sskipy = n/nlabels
	#---plot
	im = ax.imshow(array(dat).T,interpolation='nearest',origin='lower',
		norm=lognorm,cmap=cmap,alpha=cmap_washout,vmin=lims[0],vmax=lims[1])
	#---set axes labels
	if tickshow != False:
		if tickshow[0]:
			ax.set_xticks(array(list(arange(0,m/2,sskipx)*-1)[:0:-1]+list(arange(0,m/2,sskipx)))+cm)
			ax.axes.set_xticklabels([int(i) for i in array(list(arange(0,m/2,sskipx)*-1)[:0:-1]+
				list(arange(0,m/2,sskipx)))*mset.rounder/mset.lenscale])
			if ticklabel_show[0]:
				ax.set_xlabel(xlabel)
		else:
			ax.set_xticklabels([])
		if tickshow[1]:
			ax.set_yticks(array(list(arange(0,n/2,sskipy)*-1)[:0:-1]+list(arange(0,n/2,sskipy)))+cn)
			ax.axes.set_yticklabels([int(i) for i in array(list(arange(0,n/2,sskipy)*-1)[:0:-1]+
				list(arange(0,n/2,sskipy)))*mset.rounder/mset.lenscale])
			if ticklabel_show[1]:
				ax.set_ylabel(ylabel)
		else:
			ax.set_yticklabels([])
	elif tickshow == False:
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
	if inset:
		axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="30%",height="30%",loc=1)
		plotter2d(axins,mset,
			dat=data[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1],
			tickshow=False,lognorm=False,cmap=cmap,lims=lims,inset=False)
	return im

def plotter_undulate(mset,qmagfilter=None,inset2d=True,inset2d2=True,ax=None):
	'''Standard function for plotting 1D undulation spectra, with inset options.'''
	if qmagfilter == None: qmagfilter = mset.undulate_qmagfilter
	spec1d = array([i for i in array(mset.undulate_spec1d) if i[0] != 0.])
	#---some repeat comands from mset.calculate_undulations here, necessary for plotting
	m,n = shape(mset.undulate_hqhq2d)
	cm,cn = [int(i/2) for i in shape(mset.undulate_hqhq2d)]
	specfilter = array(filter(lambda x: x[0] >= qmagfilter[0] 
		and x[0] <= qmagfilter[1],mset.undulate_spec1d))
	[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
	area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] for i in mset.surf_index])/mset.lenscale**2)
	print 'kappa = '+str(1/exp(az)/area)+' kBT'
	#---calculate kappa assuming correct q4 scaling
	leftcom = [mean(log(specfilter[:,0])),mean(log(specfilter[:,1]))]
	az_enforced = leftcom[1]+4.*leftcom[0]
	kappa_enforced = 1./exp(az_enforced)/area
	print 'kappa_enforced = '+str(kappa_enforced)
	#---plot
	if ax != None:
		fig = plt.figure(figsize=(6,6))
		gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
		ax = plt.subplot(gs[0])
	ax.scatter(spec1d[:,0],spec1d[:,1],marker='o',c='k',s=20)
	ax.plot(arange(spec1d[:,0].min()/2,spec1d[:,0].max()*4),[exp(az)*(i**bz) 
		for i in arange(spec1d[:,0].min()/2,spec1d[:,0].max()*4)],linestyle='dotted',c='b',lw=2)
	ax.plot(specfilter[:,0],[exp(az)*(i**bz) for i in specfilter[:,0]],c='b',lw=2)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((spec1d[:,0].min()/2,spec1d[:,0].max()*4))
	ax.set_ylim((spec1d[:,1].min()/10,spec1d[:,1].max()*10))
	ax.set_ylabel(r'$\left\langle h_{q}h_{-q}\right\rangle$',fontsize=fsaxlabel)
	ax.set_xlabel(r'${\left|q_x\right|}(\mathrm{{nm}^{-1}})$',fontsize=fsaxlabel)
	ax.grid(True,which='major')
	labels = [item.get_text() for item in ax.get_xticklabels()]
	plt.tick_params(labelsize=fsaxticks)
	#---inset
	if inset2d:
		axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="45%",height="45%",loc=1)
		plotter2d(axins,mset)
		if inset2d2:
			i2wid = 3
			axins2 = mpl_toolkits.axes_grid.inset_locator.inset_axes(axins,width="30%",height="30%",loc=1)
			plotter2d(axins2,mset,
				dat=mset.undulate_hqhq2d[cm-i2wid:cm+i2wid+1,cn-i2wid:cn+i2wid+1],
				silence=True,cmap=mpl.cm.jet)
	#---report kappa
	ax.text(0.05,0.05,r'$\boldsymbol{\kappa_{'+str('%1.1f'%bz)+'}} = '+
		str('%3.1f'%mset.undulate_kappa)+'\:k_BT$'+
		'\n'+r'$\boldsymbol{\kappa_{-4.0}} = '+str('%3.1f'%kappa_enforced)+'\:k_BT$',
		transform=ax.transAxes,fontsize=fsaxtext)
	if ax != None:
		plt.show()

def plotmov(dat,basename,figsize=None,keep_snapshots=True,cmap=None,lowres=False,smooth=None,
	xedges=None,yedges=None):
	'''Generic wrapper for making movies. Animates a 2D plot from a numpy array.'''
	import scipy.ndimage
	#---names
	sysname = 'test'
	plot_video_filebase = 'tmpdir'
	if not os.path.isdir(pickles+'/figs-'+basename):
		os.mkdir(pickles+'/figs-'+basename)
	plot_video_filebase = '/figs-'+basename+'/'+'snap-'+basename
	#---settings
	#---Nb still needs labels and axes
	if figsize == None: figsize = (5,5)
	if cmap == None: cmap = mpl.cm.jet
	if type(dat) != 'numpy.ndarray':
		dat = array(dat)
	vmax = dat.max()
	vmin = dat.min()
	print [vmin,vmax]
	#---loop over frames
	framenums = range(len(dat))
	for fr in framenums:
		print 'Plotting frame '+str(fr)
		fig = plt.figure(figsize=figsize)
		ax = plt.subplot(111)
		if smooth != None:
			ax.imshow((scipy.ndimage.filters.gaussian_filter(dat[fr],smooth)).T,
				interpolation='nearest',origin='lower',cmap=mpl.cm.jet,vmax=vmax,vmin=vmin)
		else:
			ax.imshow((dat[fr]).T,interpolation='nearest',origin='lower',
				cmap=mpl.cm.jet,vmax=vmax,vmin=vmin)
		ax.set_xticks(range(len(xedges)))
		ax.set_yticks(range(len(yedges)))
		ax.set_xticklabels(xedges)
		ax.set_yticklabels(yedges)
		plt.savefig(pickles+plot_video_filebase+'.fr.'+str('%04d'%framenums.index(fr))+'.png',
			dpi=200,bbox_inches='tight')
		plt.close(fig)
	if lowres:
		subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png',
			'-vcodec','mpeg2video','-filter:v','setpts=2.0*PTS',
			pickles+'/'+plot_video_filebase+'.lowres.mpeg'])
	else:
		subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png',
			'-vcodec','mpeg2video','-qscale','0','-filter:v','setpts=2.0*PTS',
			pickles+'/'+plot_video_filebase+'.mpeg'])
	if keep_snapshots == False:
		os.popen('rm -r -f '+pickles+'/figs-'+sysname+'-dimple-view')







#!/usr/bin/python -i

from membrainrunner import *

import subprocess
import os

#---Note that this includes tranpose corrections for the mset.surf data.

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

location = ''
execfile('locations.py')
execfile('plotter.py')

#---filenames
struct_pickles = ['pkl.structures.membrane-v700.md.part0006.360000-460000-200.pkl',
	'pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',
	'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl']
startpickle = struct_pickles[-1]
sysname = startpickle[24:-4]

#---settings
do_average_plot = False
do_video = True

#---maps labels
mapslabels = [r'$\textbf{{EXO70}\ensuremath{\times}2{\small (parallel)}}$',
	r'$\textbf{{EXO70}\ensuremath{\times}2{\small (antiparallel)}}$']
	
#---plot settings
plt.rc('font', family='sans-serif')
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
#mpl.rcParams.update({'font.style':'sans-serif'})
#mpl.rcParams.update({'font.size': 16})
#---these broke somehow, I think
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def curvcalc(z,lenscale):
	'''Calculate mean and Gaussian curvature directly.'''
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

#---load
mset = unpickle(pickles+startpickle)

#---plot average mean and Gaussian curvature
if do_average_plot:
	#---calculate mean curvature
	mset.calculate_average_surface()
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	cdat = lenscale*curvcalc(mset.surf[0],lenscale)[0]
	curvsm = []
	for i in range(len(mset.surf)):
		curvsm.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[0])
	#---calculate Gaussian curvature	
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	curvsk = []
	for i in range(len(mset.surf)):
		curvsk.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[1])
	#---plots
	extremum = max([np.max(curvsm),np.abs(np.min(curvsm))])
	extremum = abs(np.mean(curvsm))+2*np.std(curvsm)
	numgridpts = shape(curvsk)[1]
	fig = plt.figure(figsize=(12,3))
	ax0 = plt.subplot2grid((1,3), (0,0))
	ax0.set_title('mean curvature')
	img0 = ax0.imshow(array(np.mean(curvsm,axis=0)).T,interpolation='nearest',origin='LowerLeft',
		vmax=extremum,vmin=-extremum,
		cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	vecs = np.mean(mset.vecs,axis=0)
	for pt in np.mean(mset.protein,axis=0):
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
			alpha=1.0)
		ax0.add_patch(circ)	
	cax = inset_axes(ax0,
         width="5%",
         height="100%",
         bbox_transform=ax0.transAxes,
         bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
         loc= 1)
	fig.colorbar(img0,cax=cax)
	cax.tick_params(labelsize=8) 
	cax.set_ylabel(r'$\mathsf{H(nm^{-1})}$',fontsize=10)
	extremum = max([np.max(curvsk),np.abs(np.min(curvsk))])
	extremum = abs(np.mean(curvsk))+2*np.std(curvsk)
	ax1 = plt.subplot2grid((1,3), (0,1))
	ax1.set_title('Gaussian curvature')
	img = ax1.imshow(array(np.mean(curvsk,axis=0)).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	vecs = np.mean(mset.vecs,axis=0)
	for pt in np.mean(mset.protein,axis=0):
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
			alpha=1.0)
		ax1.add_patch(circ)	
	cax = inset_axes(ax1,
         width="5%",
         height="100%",
         bbox_transform=ax1.transAxes,
         bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
         loc= 1)
	fig.colorbar(img,cax=cax)
	cax.tick_params(labelsize=10) 
	cax.set_ylabel(r'$\mathsf{K(nm^{-2})}$',fontsize=10)
	extremum = max([np.max(mset.surf_mean),np.abs(np.min(mset.surf_mean))])
	ax2 = plt.subplot2grid((1,3), (0,2))
	ax2.set_title('structure')
	img2 = ax2.imshow(array(mset.surf_mean).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
	vecs = np.mean(mset.vecs,axis=0)
	for pt in np.mean(mset.protein,axis=0):
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
			alpha=1.0)
		ax2.add_patch(circ)	
	cax = inset_axes(ax2,
         width="5%",
         height="100%",
         bbox_transform=ax2.transAxes,
         bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
         loc= 1)
	fig.colorbar(img2,cax=cax)
	cax.tick_params(labelsize=10) 
	cax.set_ylabel(r'$\mathsf{z(nm)}$',fontsize=10)
	plt.tight_layout()
	plt.savefig(pickles+'fig-'+sysname+'-curvatures.mean.png',dpi=500,bbox_inches='tight')
	plt.show()
	plt.cla()
	
#---plot average mean and Gaussian curvature
if do_video:
	if not os.path.isdir(pickles+'/figs-'+sysname+'-curvature.snapshots'):
		os.mkdir(pickles+'/figs-'+sysname+'-curvature.snapshots')
	#---calculate mean curvature
	mset.calculate_average_surface()
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	cdat = lenscale*curvcalc(mset.surf[0],lenscale)[0]
	curvsm = []
	for i in range(len(mset.surf)):
		curvsm.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[0])
	#---calculate Gaussian curvature	
	lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
	curvsk = []
	for i in range(len(mset.surf)):
		curvsk.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[1])
	#---plots
	framecount = 0
	for fr in range(0,len(mset.surf),1):
		print 'rendering frame = '+str(fr)
		extremum = max([np.max(curvsm),np.abs(np.min(curvsm))])
		extremum = abs(np.mean(curvsm))+2*np.std(curvsm)
		numgridpts = shape(curvsk)[1]
		fig = plt.figure(figsize=(12,3))
		ax0 = plt.subplot2grid((1,3), (0,0))
		ax0.set_title('mean curvature')
		img0 = ax0.imshow(array(curvsm[fr]).T,interpolation='nearest',origin='LowerLeft',
			vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		for pt in mset.protein[fr]:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax0.add_patch(circ)	
		cax = inset_axes(ax0,
		     width="5%",
		     height="100%",
		     bbox_transform=ax0.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img0,cax=cax)
		cax.tick_params(labelsize=8) 
		cax.set_ylabel(r'$\mathsf{H(nm^{-1})}$',fontsize=10)
		extremum = max([np.max(curvsk),np.abs(np.min(curvsk))])
		extremum = abs(np.mean(curvsk))+2*np.std(curvsk)
		ax1 = plt.subplot2grid((1,3), (0,1))
		ax1.set_title('Gaussian curvature')
		img = ax1.imshow(array(curvsk[fr]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		for pt in mset.protein[fr]:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax1.add_patch(circ)	
		cax = inset_axes(ax1,
		     width="5%",
		     height="100%",
		     bbox_transform=ax1.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\mathsf{K(nm^{-2})}$',fontsize=10)
		extremum = max([np.max(mset.surf_mean),np.abs(np.min(mset.surf_mean))])
		ax2 = plt.subplot2grid((1,3), (0,2))
		ax2.set_title('structure')
		img2 = ax2.imshow(array(mset.surf[fr]).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		vecs = np.mean(mset.vecs,axis=0)
		for pt in mset.protein[fr]:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2.*1./2*1.5/64*numgridpts,color='k',
				alpha=1.0)
			ax2.add_patch(circ)	
		cax = inset_axes(ax2,
		     width="5%",
		     height="100%",
		     bbox_transform=ax2.transAxes,
		     bbox_to_anchor=(0.2, 0.1, 1.00, 0.95),
		     loc= 1)
		fig.colorbar(img2,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\mathsf{z(nm)}$',fontsize=10)
		plt.tight_layout()
		plt.savefig(pickles+'/figs-'+sysname+'-curvature.snapshots/fig%05d.png'%framecount,dpi=500,bbox_inches='tight')
		framecount += 1
		plt.cla()
		plt.close()
	subprocess.call(['ffmpeg','-i',pickles+'/figs-'+sysname+'-curvature.snapshots/fig%05d.png','-vcodec','mpeg2video','-qscale','0','-filter:v','setpts=2.0*PTS',pickles+'/vid-'+sysname+'-curvature.mpeg'])
	os.popen('rm -r -f '+pickles+'/figs-'+sysname+'-tilefilter.snapshots')


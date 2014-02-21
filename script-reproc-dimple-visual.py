#!/usr/bin/python -i

from membrainrunner import *

import os
from scipy.optimize import leastsq
import matplotlib as mpl
from pylab import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

skip = 1
framecount = None
location = ''
execfile('locations.py')

import matplotlib.gridspec as gridspec
which_brewer_colors = [0,1,2,3,4,5,6,7]
clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]
clrs2 = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams.update({'font.size': 14})

from scipy import ndimage

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

analysis_descriptors = [
	('pkl.dimple.v614-stress.md.part0002.rerun.pkl',None,
		r'$\textbf{{ENTH}\ensuremath{\times}4}$'),
	('pkl.dimple.v612-stress.md.part0003.pkl',None,
		r'$\textbf{{ENTH}\ensuremath{\times}1}$'),
	('pkl.dimple.v612-stress.md.part0003.prot-v614.pkl',None,
		r'$\textbf{{ENTH}\ensuremath{\times}1(4x area)}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.pkl',None,
		r'$\textbf{control}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v614.invert.pkl',None,
		r'$\textbf{control, invert}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v614.pkl',None,
		r'$\textbf{control shift}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v614.pkl',None,
		r'$\textbf{control shift}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift01.prot-v700.pkl',None,
		r'$\textbf{control shift}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.shift10.prot-v700.pkl',None,
		r'$\textbf{control shift}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.prot-v700.pkl',None,
		r'$\textbf{control}$'),
	('pkl.dimple.v700.md.part0002.100000-200000-200.pkl',None,
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (para)}}$'),
	('pkl.dimple.v701.md.part0003.60000-160000-200.pkl',None,
		r'$\textbf{{EXO70}\ensuremath{\times}2{\small (anti)}}$'),
	('pkl.dimple.v550.md.part0006.300000-400000-200.dummytest.pkl',None,
		r'$\textbf{control, invert}$')]
analyses = [analysis_descriptors[i] for i in [-1]]

show_colorbar = False
keep_snapshots = True
protein_subset_slice = slice(None,None)
fs = 10
figsize = (20,6.45)

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	z0 = 0
	#return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over desired analyses
for ad in analyses:

	#---load
	resdat = pickle.load(open(pickles+ad[0],'r'))
	startpickle = resdat[0].getnote('startpickle')
	mset = unpickle(pickles+startpickle)	
	protein_pickle = ad[1]
	sysname = startpickle[24:-4]
	print 'loaded '+startpickle

	#---find protein points
	if protein_pickle == None:
		proteins_all = array(mset.protein)
		proteins = proteins_all[:,protein_subset_slice]
	else:
		mset_protein = unpickle(pickles+protein_pickle)
		proteins_all = array(mset_protein.protein)
		proteins = proteins_all[:,protein_subset_slice]

	#---define box vectors for grid-to-real distance conversions
	griddims = mset.griddims
	cutoff_distance = resdat[0].getnote('cutoff_distance')
	zfilterdir = resdat[0].getnote('expected_direction')
	params = []
	for zfilterdir in [-1,1,0]:
		rd = [-1,1,0].index(zfilterdir)
		params.append(resdat[rd].get(['type','params']))
	
	#---prepare to render
	if not os.path.isdir(pickles+'/figs-'+sysname+'-dimple-view'):
		os.mkdir(pickles+'/figs-'+sysname+'-dimple-view')
	plot_video_filebase = '/figs-'+sysname+'-dimple-view/'+'snap-dimple-view'

	#---frame selection header
	end = None
	start = None
	nframes = len(mset.surf)
	if framecount == None:
		if end == None: end = nframes
		if start == None: start = 0
		if skip == None: skip = 1
	else:
		start = 0
		end = nframes
		skip = int(float(nframes)/framecount)
		skip = 1 if skip < 1 else skip
	print 'frame count = '+str(end)
	framenums = range(start,end,skip)
	
	#---loop over frames
	for fr in range(start,end,skip):
	
		#---prepare figure
		fig = plt.figure(figsize=figsize)
		gs = gridspec.GridSpec(2,6,wspace=0.0,hspace=0.0)
	
		#---loop over possible height filters
		for zfilterdir in [-1,1,0]:
	
			print 'Fitting frame '+str(fr)
			vecs = mset.vec(fr)
			cutoff = cutoff_distance*10/(vecs[0]/griddims[0])
			rd = [-1,1,0].index(zfilterdir)
		
			#---MAJOR CODE BLOCK FOLLOWS

			#---tesselate (via PBCs) the original surface
			#---surfpbc = numpy.tile(mset.surf[fr],(3,3))
			surfpbc = mset.surf[fr]
			#---positive/negative domain selection
			surfpoz = array([[(1 if surfpbc[i][j] > 0 else 0) 
				for j in range(1*(griddims[1]-1))] for i in range(1*(griddims[0]-1))]).T
			surfneg = array([[(1 if surfpbc[i][j] < 0 else 0) 
				for j in range(1*(griddims[1]-1))] for i in range(1*(griddims[0]-1))]).T
			label_im_poz, nb_labels = ndimage.label(surfpoz)
			label_im_neg, nb_labels = ndimage.label(surfneg)
			surf_discrete = label_im_poz - label_im_neg
			#---protein shadow selection
			prot_disc_pts = [(int(round(pt[0]/vecs[0]*(griddims[0]-1))),int(round(pt[1]/vecs[1]*(griddims[1]-1)))) 
				for pt in mset.protein[fr]]
			#---the following is equivalent
			#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
			prot_locs = [[(1 if (i,j) in prot_disc_pts else 0) for i in range(griddims[0]-1)] 
				for j in range(griddims[1]-1)]
			#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
			gridpos = [[i,j] for i in range(griddims[0]-1) for j in range(griddims[1]-1)]
			distmat = scipy.spatial.distance.cdist(prot_disc_pts,gridpos)
			bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
				for i in range(len(distmat))] for w in v]))]
			shadow = [[(1 if [i,j] in bufferlist else 0) for i in range(griddims[0]-1)] 
				for j in range(griddims[1]-1)]
			#---final selection
			targetbase = [list(i) 
				for i in array(where((array(surf_discrete).T>0.)*(array(shadow).T==1)!=0)).T]
			if zfilterdir == 1:
				supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in targetbase])) if i > 0]
			elif zfilterdir == -1:
				supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in targetbase])) 
					if i < 0]
			else:
				supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in targetbase]))]
			if zfilterdir == 0:
				fintarget = array(bufferlist)
			else:
				fintarget = array([list([l[1],l[0]]) 
					for l in [j for k in [list(array(where(surf_discrete==i)).T) 
				for i in supplement_domains] for j in k]])
			#---convert to xyz
			targetxyz = array([[i[0]*vecs[0]/(griddims[0]-1),i[1]*vecs[1]/(griddims[1]-1),
				mset.surf[fr][i[0],i[1]]] for i in fintarget])
	
			#---calculate fitted height map
			zmap = array([[gauss2d(params[rd][fr],x*vecs[0]/(griddims[0]-1),y*vecs[1]/(griddims[1]-1)) 
				for y in range(griddims[1]-1)] for x in range(griddims[0]-1)]).T
			maxzpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
				for j in [targetxyz[argmax([gauss2d(params[rd][fr],i[0],i[1]) for i in targetxyz])]]][0]
			minzpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
				for j in [targetxyz[argmin([gauss2d(params[rd][fr],i[0],i[1]) for i in targetxyz])]]][0]
			minz = np.min([0.1*gauss2d(params[rd][fr],i[0],i[1]) for i in targetxyz])
			maxz = np.max([0.1*gauss2d(params[rd][fr],i[0],i[1]) for i in targetxyz])

			#---calculate curvature map
			curvmap = array([[-1*gauss2dh(params[rd][fr],x*vecs[0]/(griddims[0]-1),y*vecs[1]/(griddims[1]-1)) 
				for y in range(griddims[1]-1)] for x in range(griddims[0]-1)]).T
			maxhpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
				for j in [targetxyz[argmax([-1*gauss2dh(params[rd][fr],i[0],i[1]) for i in targetxyz])]]][0]
			minhpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
				for j in [targetxyz[argmin([-1*gauss2dh(params[rd][fr],i[0],i[1]) for i in targetxyz])]]][0]
			minh = np.min([-10*gauss2dh(params[rd][fr],i[0],i[1]) for i in targetxyz])
			maxh = np.max([-10*gauss2dh(params[rd][fr],i[0],i[1]) for i in targetxyz])

			#---calculate residual map
			residmap = zmap-mset.surf[fr][:-1,:-1]
			maxresid = 0.1*np.max(residmap)
			minresid = 0.1*np.min(residmap)

			domaincountlim = max(np.max(surf_discrete),abs(np.min(surf_discrete)))
			ax = fig.add_subplot(gs[0,0+2*[-1,1,0].index(zfilterdir)])
			img = ax.imshow(surf_discrete,extent=None,origin='LowerLeft',
				interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-domaincountlim,vmax=domaincountlim)
			for pt in fintarget:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
			for pt in bufferlist:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./5/100*griddims[0],color='w',alpha=1.))
			for pt in [list(i) for i in array(where(array(prot_locs)==1)).T]:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./3*2/100*griddims[0],color='k',alpha=1.))
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_title('DOMAINS',fontsize=fs)

			ax = fig.add_subplot(gs[1,0+2*[-1,1,0].index(zfilterdir)])
			curvmax = np.max(abs(curvmap))
			img = ax.imshow(curvmap,extent=None,origin='LowerLeft',
				interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-curvmax,vmax=curvmax)
			for pt in fintarget:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
			ax.add_patch(plt.Circle(maxhpos,radius=2./100*griddims[0],color='k',alpha=1.))
			ax.add_patch(plt.Circle(minhpos,radius=2./100*griddims[0],color='k',alpha=1.))
			ax.add_patch(plt.Circle(maxhpos,radius=1./100*griddims[0],color='r',alpha=1.))
			ax.add_patch(plt.Circle(minhpos,radius=1./100*griddims[0],color='b',alpha=1.))
			if show_colorbar:
				divider = make_axes_locatable(ax)
				cax = divider.append_axes("right", size="5%", pad=0.05)
				plt.colorbar(img,cax=cax)
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xlabel(r'$\mathrm{H\,({nm}^{-1})\in['+str('%3.3f'%minh)+','+str('%3.3f'%maxh)+']}$',fontsize=fs)

			ax = fig.add_subplot(gs[1,1+2*[-1,1,0].index(zfilterdir)])
			zmax = np.max(abs(zmap))
			img = ax.imshow(zmap,extent=None,origin='LowerLeft',
				interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-zmax,vmax=zmax)
			for pt in fintarget:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
			ax.add_patch(plt.Circle(maxzpos,radius=2./100*griddims[0],color='k',alpha=1.))
			ax.add_patch(plt.Circle(minzpos,radius=2./100*griddims[0],color='k',alpha=1.))
			ax.add_patch(plt.Circle(maxzpos,radius=1./100*griddims[0],color='r',alpha=1.))
			ax.add_patch(plt.Circle(minzpos,radius=1./100*griddims[0],color='b',alpha=1.))
			if show_colorbar:
				divider = make_axes_locatable(ax)
				cax = divider.append_axes("right", size="5%", pad=0.05)
				plt.colorbar(img,cax=cax)
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xlabel(r'$\mathrm{z\,(nm)\in['+str('%3.3f'%minz)+','+str('%3.3f'%maxz)+']}$',fontsize=fs)

			ax = fig.add_subplot(gs[0,1+2*[-1,1,0].index(zfilterdir)])
			residmax = np.max(abs(residmap))
			img = ax.imshow(residmap,extent=None,origin='LowerLeft',
				interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-residmax,vmax=residmax)
			for pt in fintarget:
				ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
			if show_colorbar:
				divider = make_axes_locatable(ax)
				cax = divider.append_axes("right", size="5%", pad=0.05)
				plt.colorbar(img,cax=cax)
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_title(r'$\mathrm{(z-\hat{z})\,(nm)\in['+str('%3.1f'%minresid)+','+str('%3.1f'%maxresid)+']}$',fontsize=fs)
		plt.savefig(pickles+plot_video_filebase+'.fr.'+str('%04d'%framenums.index(fr))+'.png',dpi=200,bbox_inches='tight')
		plt.cla()
		plt.close()
	subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png','-vcodec','mpeg2video','-qscale','0','-filter:v','setpts=2.0*PTS',pickles+'/'+plot_video_filebase+'.mpeg'])
	subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png','-vcodec','mpeg2video','-filter:v','setpts=2.0*PTS',pickles+'/'+plot_video_filebase+'.lowres.mpeg'])
	if keep_snapshots == False:
		os.popen('rm -r -f '+pickles+'/figs-'+sysname+'-dimple-view')


#!/usr/bin/python -i

from membrainrunner import *

import os
from scipy.optimize import leastsq
import matplotlib as mpl
from pylab import *
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---parameters
skip = 1
framecount = None
location = ''
execfile('locations.py')

#---analysis plan
analysis_plan = slice(-6,None)
analysis_descriptors = [
	('pkl.structures.membrane-v701.md.part0003.60000-160000-200.pkl',slice(None),None,-1,False,''),
	('pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',slice(None),None,-1,False,''),
	('pkl.structures.membrane-v612-stress.md.part0003.pkl',slice(None),None,1,False,''),
	('pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',slice(None),None,1,False,''),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.prot-v614'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,(0,1),'.shift01.prot-v614'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,(1,0),'.shift10.prot-v614'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',1,False,'.prot-v700'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',1,(0,1),'.shift01.prot-v700'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v700.md.part0002.100000-200000-200.pkl',1,(1,0),'.shift10.prot-v700'),
	('pkl.structures.membrane-v612-stress.md.part0003.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.prot-v614'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.prot-v614.invert'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.dummytest'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.prot-v614.phasetest'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[16,16],'.prot-v614.shift-16-16'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[8,8],'.prot-v614.shift-8-8'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[24,24],'.prot-v614.shift-24-24'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[16,0],'.prot-v614.shift-16-0'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[0,16],'.prot-v614.shift-0-16'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[1,1],'.prot-v614.shift-1-1'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[3,3],'.prot-v614.shift-3-3'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[4,4],'.prot-v614.shift-4-4'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[5,5],'.prot-v614.shift-5-5'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[2,0],'.prot-v614.shift-2-0'),
	('pkl.structures.membrane-v550.md.part0006.300000-400000-200.pkl',slice(None),
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,[0,2],'.prot-v614.shift-0-2')]
	
'''
NOTES NOTES NOTES

other parameters:
XY filter 	: cutoff_distance 
Z filter 	: positive or negative or none
contiguity	: keep everything with a positive
Hmax filter	: which range of Hmax
sig filter	: which extents

whether to look just at the Hmaxs or the whole distribution?
whether to include dimples with centers in the right region, or anywhere?
if we don't want points in the shadow, then can we constrain the optimizer?
what position, size, shape should the shadow take when we study the control?

current plan:

'''

#---parameters
cutoff_distance = 15.

#---methods
special_inversion_test = False
special_dummy_test = False
special_phase_test = False
special_shift_test = True

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---Load expressions for mean and Gaussian curvature
execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	#---fix the height shift
	z0 = 0
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2d_residual(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual.'''
	return ((gauss2d(params,x,y)-z)**2)
	
def lateral_discretize(fr,result='area'):
	'''
	LATERAL_DISCRETIZE
	Discretize a boolean of bilayer midplane sign given access to a membrane object and a proteins object.
	Return protein present/absent region.
	Options:
	result : default return just the area calculation, grid returns the result grids, all returns both
	'''
	print 'Analyzing discretized areas for frame: '+str(fr)
	#---transpose correction
	surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) for j in range(mset.griddims[1]-1)] 
		for i in range(mset.griddims[0]-1)]).T
	protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(protpos,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
		dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
	if result == 'all':
		return [[float(sum(surf_discrete==1)),float(sum(surf_discrete==-1)),
			float(sum(surf_discrete+buf==2)),float(sum(surf_discrete+buf==0)),float(sum(buf==1))],
			surf_discrete,buf]
	elif result == 'grid':
		return [surf_discrete,buf]
	else:
		return [float(sum(surf_discrete==1)),float(sum(surf_discrete==-1)),
			float(sum(surf_discrete+buf==2)),float(sum(surf_discrete+buf==0)),float(sum(buf==1))]
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over desired analyses
for ad in analysis_descriptors[analysis_plan]:
	(startpickle,protein_subset_slice,protein_pickle,expected_direction,testshift,suffix) = ad
	sysname = startpickle[24:-4]

	#---load
	mset = unpickle(pickles+startpickle)
	print 'loaded '+startpickle
	result_data_collection = []

	#---find protein points
	if protein_pickle == None:
		proteins_all = array(mset.protein)
		proteins = proteins_all[:,protein_subset_slice]
	else:
		mset_protein = unpickle(pickles+protein_pickle)
		proteins_all = array(mset_protein.protein)
		proteins = proteins_all[:,protein_subset_slice]

	#---define box vectors for grid-to-real distance conversions
	griddims = mset.griddims[0]
	vecs = np.mean(mset.vecs,axis=0)
	cutoff = cutoff_distance*10/(vecs[0]/mset.griddims[0])

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
	
	#---loop over possible surface point selection filters
	for zfilterdir in [-1,1,0]:
		#---fit and save

		#---DEPRECATED
		#[params,maxhs,maxhxys,target_zones,which_frames] = batch_dimple_fitting(skip=None,
		#	framecount=framecount)
		
		#---should you fix the vecs so it's framewise? and change mset.griddims to griddims for clarity

		#---declare
		params = []
		maxhs = []
		maxhxys = []
		target_zones = []
		residsum = []
		
		#---loop over frames
		for fr in range(start,end,skip):
			print 'Fitting frame '+str(fr)
		
			#---MAJOR CODE BLOCK FOLLOWS

			#---frame-specific cutoff
			griddims = mset.griddims
			#---tesselate (via PBCs) the original surface
			#---surfpbc = numpy.tile(mset.surf[fr],(3,3))
			if special_inversion_test:
				surfpbc = list(-1*array(mset.surf[fr]))
			elif special_phase_test:
				randshift = [randint(0,31) for i in range(2)]
				surfpbc = array([[mset.surf[fr][(j+randshift[0])%(mset.griddims[0]-2)]\
					[(i+randshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))]).T
			elif special_shift_test:
				surfpbc = array([[mset.surf[fr][(j+testshift[0])%(mset.griddims[0]-2)]
					[(i+testshift[1])%(mset.griddims[1]-2)] for j in range(1*(mset.griddims[1]-1))] 
					for i in range(1*(mset.griddims[0]-1))]).T
			else:
				surfpbc = mset.surf[fr]
			#---positive/negative domain selection
			surfpoz = array([[(1 if surfpbc[i][j] > 0 else 0) 
				for j in range(1*(mset.griddims[1]-1))] for i in range(1*(mset.griddims[0]-1))]).T
			surfneg = array([[(1 if surfpbc[i][j] < 0 else 0) 
				for j in range(1*(mset.griddims[1]-1))] for i in range(1*(mset.griddims[0]-1))]).T
			label_im_poz, nb_labels = ndimage.label(surfpoz)
			label_im_neg, nb_labels = ndimage.label(surfneg)
			surf_discrete = label_im_poz - label_im_neg
			#---protein shadow selection
			prot_disc_pts = [[int(round(pt[0]/vecs[0]*(griddims[0]-1))),g
				int(round(pt[1]/vecs[1]*(griddims[1]-1)))] for pt in proteins[fr%len(proteins)]]
			#---the following is equivalent
			#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
			prot_locs = [[(1 if [i,j] in prot_disc_pts else 0) for i in range(mset.griddims[0]-1)] 
				for j in range(mset.griddims[1]-1)]
			#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
			gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
			distmat = scipy.spatial.distance.cdist(prot_disc_pts,gridpos)
			bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
				for i in range(len(distmat))] for w in v]))]
			shadow = [[(1 if [i,j] in bufferlist else 0) for i in range(mset.griddims[0]-1)] 
				for j in range(mset.griddims[1]-1)]
			#---final selection
			targetbase = [list(i) 
				for i in array(where((array(surf_discrete).T>0.)*(array(shadow).T==1)!=0)).T]
			if zfilterdir == 1:
				supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in targetbase])) 
					if i > 0]
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
			if len(fintarget) > 0:
				#---convert to xyz
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][i[0],i[1]]] for i in fintarget])
				target_com = [mean(fintarget[:,0]),mean(fintarget[:,1])]

				#---Perform the fit
				if len(target) >= 7:
					p_opt = leastsq(gauss2d_residual,array([0,1,target_com[0],target_com[0],50,50,0]),
						args=(target[:,0],target[:,1],target[:,2]))
					params.append(p_opt[0])

					#---find the maximum curvature strength position
					maxhxys.append(argmax([abs(gauss2dh(p_opt[0],i[0],i[1])) for i in target]))

					#---save the residuals
					residsum.append(sum([abs(gauss2d(p_opt[0],i[0],i[1])) for i in target]))

					#---reverse curvature sign here
					maxhs.append(-1*gauss2dh(p_opt[0],target[maxhxys[-1]][0],target[maxhxys[-1]][1]))

					#---save the target points
					target_zones.append(target)
			
					#---save the analyzed frames
					which_frames = range(start,end,skip)

		#---save the data to the MembraneData object			
		result_data = MembraneData('dimple',label=sysname)
		for i in range(len(params)):
			result_data.add([params[i],maxhs[i],maxhxys[i],target_zones[i]],[which_frames[i]])
		result_data.addnote(['startpickle',startpickle])
		result_data.addnote(['zfilterdir',zfilterdir])
		result_data.addnote(['cutoff_distance',cutoff_distance])
		result_data.addnote(['protein_subset_slice',protein_subset_slice])
		result_data.addnote(['protein_pickle',protein_pickle])
		result_data.addnote(['expected_direction',expected_direction])
		result_data.addnote(['sysname',sysname])
		result_data.addnote(['suffix',suffix])
		result_data_collection.append(result_data)
		del result_data

	#---write the data
	pickle.dump(result_data_collection,open(pickles+'pkl.dimple.'+sysname+suffix+'.pkl','w'))		

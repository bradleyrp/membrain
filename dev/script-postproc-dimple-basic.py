#!/usr/bin/python -i

from membrainrunner import *

import os
from scipy.optimize import leastsq
import matplotlib as mpl
from pylab import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---parameters
skip = 1
framecount = None
location = ''
execfile('locations.py')

#---analysis plan
analysis_plan = slice(-2,-1)
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
		'pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',1,False,'.dummytest'),]
	
#---parameters
cutoff_distance = 15.
curvature_filter = [0.001,0.1]
special_inversion_test = True
special_dummy_test = False

#---plot settings
original_plot_style = False
which_brewer_colors = [0,2,3]
colorcodes = [brewer2mpl.get_map('Set2','qualitative',8).mpl_colors[i] for i in which_brewer_colors]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---Load expressions for mean and Gaussian curvature
execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	z0 = 0
	#return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2d_residual(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual.'''
	#return ((g2d(params,x,y)-z)**2)*(1 if center == None else 1./norm([x-center[0],y-center[1]]))
	#return ((g2d(params,x,y)-z)**2)/norm([x-center[0],y-center[1]])**2
	return ((gauss2d(params,x,y)-z)**2)
	
def batch_dimple_fitting(end=None,start=None,skip=None,framecount=None):
	'''Loop over frames, and fit height profiles.'''
	params = []
	maxhs = []
	maxhxys = []
	target_zones = []
	residsum = []
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
	cutoff = cutoff_distance*10/(vecs[0]/mset.griddims[0])
	print 'Total frames = '+str(end)
	if special_dummy_test:
		dummysurf = array([[i*vecs[0]/(mset.griddims[0]-1),j*vecs[1]/(mset.griddims[1]-1),
			gauss2d([20,10.,vecs[0]/2.,vecs[1]/2.,10,20,0],i*vecs[0]/(mset.griddims[0]-1),
			j*vecs[1]/(mset.griddims[1]-1))] for j in range(mset.griddims[1]-1) 
			for i in range(mset.griddims[0]-1)])
	for fr in range(start,end,skip):
		print 'Fitting frame '+str(fr)
		#---note: fixed midplaner transpose error here
		#---note: many transposes are just for the "where" command
		if testshift != False:
			grids = mset.griddims
			shift = [testshift[0]*int(grids[0]/2.),testshift[1]*int(grids[1]/2.)]
			surf_discrete = array([[(1 if mset.surf[fr][(i+shift[0])%grids[0]][(j+shift[1])%grids[1]] > 0 
				else -1) for j in range(mset.griddims[1]-1)] for i in range(mset.griddims[0]-1)]).T
		elif special_inversion_test == True:
			surf_discrete = array([[(1 if -mset.surf[fr][i][j] > 0 else -1) 
				for j in range(mset.griddims[1]-1)] 
				for i in range(mset.griddims[0]-1)]).T
		elif special_dummy_test == True:
			surf_discrete = array([[1 for j in range(mset.griddims[1]-1)] 
				for i in range(mset.griddims[0]-1)]).T			
		else:
			surf_discrete = array([[(1 if mset.surf[fr][i][j] > 0 else -1) 
				for j in range(mset.griddims[1]-1)] 
				for i in range(mset.griddims[0]-1)]).T
		if testshift != False:
			protpos = array(where(array(mset.rezipgrid(mean(proteins,axis=0)))!=0.0)).T
		else:
			protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
		gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
		#---find distances between protein points and all grid points.
		distmat = scipy.spatial.distance.cdist(protpos,gridpos)
		#---find grid points within cutoff of protein points
		bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
			for i in range(len(distmat))] for w in v]))]
		#---find absolute coordinates for grid points within the protein "buffer"
		buf = array(scipy.sparse.coo_matrix(([1 for i in range(len(bufferlist))],array(bufferlist).T),
			dtype=int8,shape=(mset.griddims[0]-1,mset.griddims[1]-1)).todense())
		if testshift != False:
			#---select target for fitting
			if height_direction == 1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][(i[0]+shift[0])%grids[0],(i[1]+shift[1])%grids[1]]] 
					for i in array(where(surf_discrete+buf==2)).T])
			elif height_direction == -1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][(i[0]+shift[0])%grids[0],(i[1]+shift[1])%grids[1]]] 
					for i in array(where(surf_discrete+buf==0)).T])
			elif height_direction == 0:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][(i[0]+shift[0])%grids[0],(i[1]+shift[1])%grids[1]]] 
					for i in array(where(buf==1)).T])
		elif special_dummy_test == True:
				target = dummysurf
		elif special_inversion_test == True:
			#---select target for fitting
			if height_direction == 1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					-mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete+buf==2)).T])
			elif height_direction == -1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					-mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete+buf==0)).T])
			elif height_direction == 0:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					-mset.surf[fr][i[0],i[1]]] for i in array(where(buf==1)).T])
		else:
			#---select target for fitting
			if height_direction == 1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete+buf==2)).T])
			elif height_direction == -1:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][i[0],i[1]]] for i in array(where(surf_discrete+buf==0)).T])
			elif height_direction == 0:
				target = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
					mset.surf[fr][i[0],i[1]]] for i in array(where(buf==1)).T])
		#---Identify center of target points for possible location weighting
		#---Note: sometimes fitting something that isn't x,y specific (like H), you can use protein itself
		target_com = [mean(target[:,0]),mean(target[:,1])]
		#---Perform the fit
		p_opt = leastsq(gauss2d_residual,array([0,1,target_com[0],target_com[0],50,50,0]),
			args=(target[:,0],target[:,1],target[:,2]))
		params.append(p_opt[0])
		#---Store the location of maximum mean curvature and target zone
		#---Find the maximum curvature strength, regardless of direction, and save it.
		maxhxys.append(argmax([abs(gauss2dh(p_opt[0],i[0],i[1])) for i in target]))
		residsum.append(sum([abs(gauss2d(p_opt[0],i[0],i[1])) for i in target]))
		#---Reverse curvature sign here (I think I did this to match the meaning of "positive curvature")
		maxhs.append(-1*gauss2dh(p_opt[0],target[maxhxys[-1]][0],target[maxhxys[-1]][1]))
		target_zones.append(target)
	return [params,maxhs,maxhxys,target_zones,residsum,range(start,end,skip)]
	
def view_figures_timeseries_original(params=None,maxhs=None,maxhxys=None,target_zones=None):
	'''Plot the framewise measurements and residuals.'''
	#---Calculate point-wise residuals
	resids = [mean([abs(gauss2d(params[0],i[0],i[1])-i[2]) for i in target_zones[j]]) 
		for j in range(len(params))]
	validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] > curvature_filter[0] 
		and 10*maxhs[i] < curvature_filter[1])]
	#---nanometer correction in the calculation
	validhs = [10*maxhs[i] for i in validhis]
	sigma_x = [abs(params[i][4])/10. for i in validhis]
	sigma_y = [abs(params[i][5])/10. for i in validhis]
	nfitpts = [len(i) for i in target_zones]
	mpl.rcParams.update({'font.size': 30})
	fig = plt.figure(figsize=(10,10))
	ax1 = plt.subplot2grid((4,1),(0,0))
	ax1.plot(range(len(validhis)),validhs)
	ax1.set_xlim((0,len(validhis)))
	ylabel1 = ax1.set_ylabel('$H_{max}$')
	ax1.grid()
	ax2 = plt.subplot2grid((4,1),(1,0))
	ax2.plot(range(len(validhis)),sigma_x)
	ax2.plot(range(len(validhis)),sigma_y)
	ax2.set_yscale('log')
	ax2.set_xlim((0,len(validhis)))
	ylabel2 = ax2.set_ylabel('($\sigma_x,\sigma_y)$')
	ax2.grid()
	ax3 = plt.subplot2grid((4,1),(2,0))
	ax3.plot(range(len(validhis)),[resids[i] for i in validhis])
	ax3.set_xlim((0,len(validhis)))
	ylabel3 = ax3.set_ylabel('$\sum(z-\hat{z})$')
	ax3.grid()
	ax4 = plt.subplot2grid((4,1),(3,0))
	ax4.plot(range(len(validhis)),[nfitpts[i] for i in validhis])
	ax4.set_xlim((0,len(validhis)))
	ax4.set_ylabel('$N_{points}$')
	ylabel4 = ax4.set_xlabel('Frame')
	ax4.grid()
	plt.savefig(pickles+'fig-'+sysname+'-dimple.framewise.png',dpi=300,
		bbox_extra_artists=[ylabel1,ylabel2,ylabel3,ylabel4],bbox_inches='tight')
	plt.close()

def view_figures_curvatures_original(params=None,maxhs=None,maxhxys=None,
	means='linear',scales='linear',fullrange=True,extrasuffix='',sigma_calc='mode'):
	'''Summarize the results in figure saved to a file. Original version.'''
	if fullrange == True:
		#---Nanometer correction
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
	else:
		#---Nanometer correction
		if height_direction == 0:
			print 'hdir 0'
			validhis = [i for i in range(len(maxhs)) if (abs(10*maxhs[i]) > curvature_filter[0] 
				and abs(10*maxhs[i]) < curvature_filter[1])]
		elif height_direction == 1:
			print 'hdir 1'
			validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] > curvature_filter[0] 
				and 10*maxhs[i] < curvature_filter[1])]
		elif height_direction == -1:
			print 'hdir -1'
			validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] < -1*curvature_filter[0] 
				and 10*maxhs[i] > -1*curvature_filter[1])]
	print validhis
	print maxhs
	#---Nanometer correction
	validhs = [10*maxhs[i] for i in validhis]
	#mpl.rcParams.update({'font.size': 30})
	#mpl.rcParams.update({'font.style':'sans-serif'})
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((2,1),(0,0))
	#---Nanometer correction in the calculation
	if scales == 'log':
		ax.hist([log10(i) for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,linewidth=2)
	elif scales == 'linear':
		ax.hist([i for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,bins=40,linewidth=2)		
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	ax.set_xlim((-0.10,0.10))
	if means == 'linear':
		mean_validhs = mean(validhs)
	elif means == 'log':
		mean_validhs = exp(mean(log(validhs)))
	if scales == 'log':
		plt.xlabel('Log $H_{max}$ $(nm^{-1})$', labelpad = 10)
	elif scales == 'linear':
		plt.xlabel(r'$H_{max}\:(\textsf{nm}^{-1})$', labelpad = 10)	
	#ax.text(0.65,0.75,r'$\left\langle H_{max}\right\rangle='+str('%.3f'%mean_validhs)+
	#	r'\:\textsf{nm}^{-1}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+'$ frames',
	#	transform=ax.transAxes,fontsize=28)
	ylabel1 =plt.ylabel('Frames', labelpad = 10)
	#plt.title('Fitted $H_{max}$, within $'+str(cutoff_distance)+'$ nm of protein')
	legend1 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax.grid()
	ax2 = plt.subplot2grid((2,1),(1,0))
	#---Nanometer correction below
	if means == 'linear':
		sigma_x = [abs(params[i][4])/10. for i in validhis]
		sigma_y = [abs(params[i][5])/10. for i in validhis]
	elif means == 'log':
		sigma_x = [log10(abs(params[i][4])/10.) for i in validhis]
		sigma_y = [log10(abs(params[i][5])/10.) for i in validhis]
	if scales == 'log':
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.8,histtype='stepfilled',linewidth=2)
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.8,histtype='stepfilled',linewidth=2)
	elif scales == 'linear':
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.8,histtype='stepfilled',bins=250,linewidth=2)
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.8,histtype='stepfilled',bins=250,linewidth=2)
		ax2.set_xlim((0,50))
		print 'Extent data are really skewed. Consider norming this one.'
		print 10**mean([log10(abs(params[i][4])/10.) for i in validhis])
		print 10**mean([log10(abs(params[i][5])/10.) for i in validhis])
	#ax2.spines['top'].set_visible(False)
	#ax2.spines['right'].set_visible(False)
	#ax2.spines['bottom'].set_position(('outward', 20))
	#ax2.spines['left'].set_position(('outward', 30))
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	for tick in ax2.xaxis.get_major_ticks():
		tick.label.set_fontsize(22)
	for tick in ax2.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	if sigma_calc == 'mean':
		sigma_a_calc = mean(sigma_x)
		sigma_b_calc = mean(sigma_y)
	elif sigma_calc == 'logmean':
		sigma_a_calc = exp(mean(sigma_x))
		sigma_b_calc = exp(mean(sigma_y))
	elif sigma_calc == 'mode':
		hist1 = histogram(sigma_x,bins=250)
		hist2 = histogram(sigma_y,bins=250)
		sigma_a_calc = hist1[1][argmax(hist1[0])]
		sigma_b_calc = hist2[1][argmax(hist2[0])]
	if scales == 'log-transformed':
		plt.xlabel('Log extent $(\sigma_a,\sigma_b)$ $\log_{10}$(nm)', labelpad = 10)
		#ax2.text(0.05,0.65,r'$\sigma_a='+str('%.3f'%sigma_a_calc)+'$ $\\textsf{nm}$'+'\n'+r'$\sigma_b='+
		#	str('%.3f'%sigma_b_calc)+'$ $\\textsf{nm}$',transform=ax2.transAxes,fontsize=14)
	elif scales == 'linear':
		plt.xlabel(r'Extent $(\sigma_a,\sigma_b)\:(\textsf{nm})$', labelpad = 10)
		#ax2.text(0.65,0.75,r'$\sigma_a='+str('%.3f'%sigma_a_calc)+r'\:\textsf{nm}$'+'\n'+r'$\sigma_b='+
		#	str('%.3f'%sigma_b_calc)+r'\:\textsf{nm}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+
		#	'$ frames',transform=ax2.transAxes,fontsize=14)
	#plt.title('Fitted extent ($\sigma_a,\sigma_b$), within $'+str(cutoff_distance)+'$ nm of protein')
	ylabel2 = plt.ylabel('Frames', labelpad = 10)
	legend2 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax2.grid()
	plt.tight_layout()
	plt.savefig(pickles+'fig-'+sysname+'-dimple.png',dpi=500,
		bbox_extra_artists=[ylabel1,ylabel2],bbox_inches='tight')
	plt.close()

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over desired analyses
for ad in analysis_descriptors[analysis_plan]:
	(startpickle,protein_subset_slice,protein_pickle,expected_direction,testshift,suffix) = ad
	sysname = startpickle[24:-4]
	#---load
	mset = unpickle(pickles+startpickle)
	print 'loaded '+startpickle
	print 'frame count = '+str(len(mset.surf[0]))
	#---find protein points
	if protein_pickle == None:
		proteins_all = array(mset.protein)
		proteins = proteins_all[:,protein_subset_slice]
	else:
		mset_protein = unpickle(pickles+protein_pickle)
		proteins_all = array(mset_protein.protein)
		proteins = proteins_all[:,protein_subset_slice]
	vecs = np.mean(mset.vecs,axis=0)
	#---note: cutoff is global
	cutoff = cutoff_distance*10/(vecs[0]/mset.griddims[0])
	result_data_collection = []
	for height_direction in [-1,1,0]:
		#---fit and save
		[params,maxhs,maxhxys,target_zones,which_frames] = batch_dimple_fitting(skip=None,
			framecount=framecount)
		result_data = MembraneData('dimple',label=sysname)
		for i in range(len(params)):
			result_data.add([params[i],maxhs[i],maxhxys[i],target_zones[i]],[which_frames[i]])
		result_data.addnote(['height_direction',height_direction])
		result_data.addnote(['cutoff',cutoff])
		result_data.addnote(['filter_low',curvature_filter[0]])
		result_data.addnote(['filter_high',curvature_filter[1]])
		result_data.addnote(['protein_subset_slice',protein_subset_slice])
		result_data.addnote(['protein_pickle',protein_pickle])
		result_data.addnote(['expected_direction',expected_direction])
		result_data.addnote(['sysname',sysname])
		result_data.addnote(['suffix',suffix])
		result_data_collection.append(result_data)
		del result_data
	pickle.dump(result_data_collection,open(pickles+'pkl.dimple.'+sysname+suffix+'.pkl','w'))
	if original_plot_style:
		#---get results from the data object
		params = result_data_collection[(0 if expected_direction == -1 else 1)].get(['type','params'])
		maxhs = result_data_collection[(0 if expected_direction == -1 else 1)].get(['type','maxhs'])
		maxhxys = result_data_collection[(0 if expected_direction == -1 else 1)].get(['type','maxhxys'])
		target_zones = result_data_collection[(0 if expected_direction == -1 else 1)].get(['type',
			'target_zones'])
		#---plot the results
		view_figures_curvatures_original(params=params,maxhs=maxhs,maxhxys=maxhxys)
		view_figures_timeseries_original(params=params,maxhs=maxhs,maxhxys=maxhxys,target_zones=target_zones)
	else:
		print 'Skipping the plotting. Please change flags or consult the reproc script.'
		

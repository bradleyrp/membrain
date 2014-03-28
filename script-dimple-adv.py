#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

from scipy.optimize import leastsq

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---define tests
testlist_std = ['full','monomers']
testlist_mono = ['full']
testlist_control = ['peak','valley']

#---plan
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4\,(v2)}$',
		'nprots':4,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'testlist':testlist_std},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1\,(v2)}$',
		'nprots':1,
		'protein_pkl':None,
		'whichframes':slice(None,None),
		'testlist':testlist_mono},
	'v550-400000-500000-160':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc',
		'label':r'$\mathrm{control\,(v2)}$',
		'nprots':1,
		'whichframes':slice(0,500),
		'protein_pkl':'v612-75000-175000-200',
		'testlist':testlist_control+[['test',[328.,279.]]]},
	'v614-40000-140000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'trajsel':'s6-sim-lonestar/md.part0002.40000-140000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4\,(v1)}$',
		'nprots':4,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'testlist':testlist_std},
	'v612-10000-80000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'s9-trestles/md.part0003.10000-80000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1\,(v1)}$',
		'nprots':1,
		'protein_pkl':None,
		'whichframes':slice(None,None),
		'testlist':testlist_mono},
	'v550-300000-400000-200':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'label':r'$\mathrm{control\,(v1)}$',
		'nprots':1,
		'whichframes':slice(0,500),
		'protein_pkl':'v612-75000-175000-200',
		'testlist':testlist_control+[['test',[328.,279.]]]}}
analysis_names = [
	'v614-120000-220000-200',
	'v612-75000-175000-200',
	'v550-400000-500000-160',
	'v614-40000-140000-200',
	'v612-10000-80000-200',
	'v550-300000-400000-200',
	][:]
routine = ['compute','plot','plot_residmaps'][1:2]
bigname = 'v614-v612-v550-combine'

#---parameters
cutoff = 320
decay_z0 = False

#---parameter set for sweeps
params = [[c,d] for c in [50,80,100,110,120,130,150,200,250] for d in [0,-1,1]]
params_plot = [[c,d] for d in [0,1] for c in [50,80,100,120,130,150]]
params_plot = [[c,d] for d in [1] for c in [50]]
params_plot = [[c,d] for c in [50,80,100,110,120,130,150,200,250] for d in [0,-1,1]]
params_plot = [[c,d] for d in [0,1] for c in [100]]
decay_z0 = True

#---plot settings
hifilt = 0.05
smallfilt = 0.001
extrem_c0_limit = 10.
hmax_nbins = 31
show_plots = False
inset_raw_hmax = False
inset_extents = True
extent_range = 16

#---under construction
residplot = False

#---filters
resid_filters = [(6,0.25),(8,0.5)]

'''
NOTES:
This script handles all fitting and plotting for the dimple-fitting algorithm.
Each simulation has a pickle containing a list of result_data objects from different fits.
Each result_data object has a notes section with the parameters for the fit.
The parameters that we may wish to sweep includes the neighborhood cutoff, and the optional z>0, z<0 filters
It is also possible to change whether the curves decay to z = 0 at long distances.
'''

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---load expressions for mean and Gaussian curvature
execfile('script-curvature-calculus.py')

def gauss2d(params,x,y,z0=True):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2d_residual(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual.'''
	return ((gauss2d(params,x,y)-z)**2)
	
def gauss2d_z0(params,x,y,z0=True):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, decay height fixed.'''
	z0,c0,x0,y0,sx,sy,th = params
	z0 = 0
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2d_residual_z0(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual, decay height fixed'''
	return ((gauss2d_z0(params,x,y)-z)**2)
	
#---hack to enable older pickles
extra_locs = [
	'./','structures-broken-transposer-error/',
	'backup-2013.20.21-enth-review-pickles/']
def unpickle_special(name):
	'''Custom un-pickle code to search for previous pickles.'''
	mset = None
	for loc in extra_locs:
		mset = unpickle(pickles+loc+name)
		if mset != None: break
	return mset
	
#---sensible colors library
def colordict(ask):
	'''Color lookup for plotting protein hulls.'''
	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in [0,1,2,3,4,6,7,8]]
	clrsdict_specific = {'red':clrs[0],'blue':clrs[1],'black':'k'}
	clrsdict_leftovers = list(array(clrs)[range(2,len(clrs))])
	if type(ask) == str:
		return clrsdict_specific[ask]
	elif type(ask) == int:
		return clrsdict_leftovers[ask%len(clrsdict_leftovers)]		
	else: return 'y'
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		anum = analysis_names.index(aname)
		mset = unpickle_special('pkl.structures.'+filespec+'.pkl')
		msets[anum] = mset

#---fits
if 'compute' in routine:
	#---load the current set of msdats which may already have some of the requested tests in the sweep
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		#---load the pickle to see if it exists	
		filespec = specname_guess(sysname,trajsel)
		pklname = 'pkl.dimple2.'+filespec+'.pkl'
		msdat = unpickle_special('pkl.dimple2.'+filespec+'.pkl')
		if msdat == None: msdats[anum] = []
		else: msdats[anum] = msdat
	#---outer loop is over parameters in the sweep
	for param in params:
		cutoff,zfiltdir = param
		#---next loop is over systems with new fits appended the msdats list and dumped to original pickle
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			filespec = specname_guess(sysname,trajsel)
			pklname = 'pkl.dimple2.'+filespec+'.pkl'
			#---load surface points			
			mset = msets[anum]
			#---protein accounting
			if protein_pkl == None: mset_protein = msets[anum]
			#---retrieve protein points from a different dictionary item
			else:
				protein_pkl_name = protein_pkl
				for i in analysis_descriptors[protein_pkl_name]: 
					vars()[i] = (analysis_descriptors[protein_pkl_name])[i]
				filespec = specname_guess(sysname,trajsel)
				mset_protein = unpickle_special('pkl.structures.'+filespec+'.pkl')
				#---re-populate the variables for the original system now that you got the protein pickle
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			nresprot = shape(mset_protein.protein)[1]/nprots
			prot_slices = [slice(i*nresprot,(i+1)*nresprot) for i in range(nprots)]
			#---innermost loop is over the following list of geographic locations relative to proteins
			tl = []
			for test in testlist:
				if test == 'full':
					tl.append(slice(None,None))
				elif test == 'monomers':
					for ps in prot_slices: tl.append(ps)
				elif test == 'valley':
					minpos = unravel_index(mean(mset.surf,axis=0).argmin(),mset.surf[0].shape)
					minxy = [minpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
					protein_com = mean(mean(mset_protein.protein,axis=0),axis=0)[:2]
					shift = minxy - protein_com
					tl.append(shift)
				elif test == 'peak':
					maxpos = unravel_index(mean(mset.surf,axis=0).argmax(),mset.surf[0].shape)
					maxxy = [maxpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
					protein_com = mean(mean(mset_protein.protein,axis=0),axis=0)[:2]
					shift = maxxy - protein_com
					tl.append(shift)
				elif type(test) == list:
					protein_com = mean(mean(mset_protein.protein,axis=0),axis=0)[:2]
					shift =  test[1] - protein_com
					tl.append(shift)
			#---if monomer-specificity is requested, replace the names
			if 'monomers' in testlist:
				ind_cleave = testlist.index('monomers')
				testlist = testlist[:ind_cleave]+\
				['monomer '+str(i) for i in range(nprots)]+\
				testlist[ind_cleave+1:]
			#---check if the tests are already in the pickle
			print 'status: running cut = '+str(cutoff)+' filter = '+str(zfiltdir)+' on '+aname	
			chopblock = []
			for k in range(len(tl)):
				test = tl[k]
				index = [i for i in range(len(msdats[anum])) if all([
				msdats[anum][i].getnote('cutoff') == cutoff,
					msdats[anum][i].getnote('zfiltdir') == zfiltdir,
					all(msdats[anum][i].getnote('tl')[msdats[anum][i].getnote('this_test')] == test)])]
				for i in index: chopblock.append(i)
			tl = [tl[i] for i in range(len(tl)) if i not in chopblock]
			if tl == []:
				print 'status: no subjects remaining for this test'
			else:
				#---perform fits over protein slices
				print 'status: running tests with tl = '+str(tl)
				params = []
				maxhs = []
				maxhxys = []
				maxhxys = []
				target_zones = []
				for fr in range(len(mset.surf))[whichframes]:
					if (fr%100) == 0: print fr
					print fr
					params_frame = []
					maxhs_frame = []
					maxhxys_frame = []
					target_zones_frame = []
					for ps in tl:
						#---define so-called protein points
						if type(ps) == slice:
							protpts = mset_protein.protein[fr][ps,0:2]
						else: 
							protpts = mset_protein.protein[fr][:,0:2]+ps
						#---get surface points
						surfpts = mset.wrappbc(mset.unzipgrid(mset.surf[fr],vecs=mset.vec(0)),
							vecs=mset.vec(fr),mode='nine')
						#---find minimum distances
						cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2])
						tmp = array([np.min(cd,axis=0),surfpts[:,2]]).T
						selected = where(tmp[:,0]<cutoff)[0]
						target = surfpts[selected]
						if zfiltdir == 1:
							target = target[target[:,2]>0]
						elif zfiltdir == -1:
							target = target[target[:,2]<0]
						#---find the center of the target as initial guess for the fitter
						target_com = [mean(target[:,0]),mean(target[:,1])]
						#---perform the fit
						#---recall z0,c0,x0,y0,sx,sy,th = params
						#---check if there are any valid points
						if len(target) >= 7:
							p_opt = leastsq((gauss2d_residual_z0 if decay_z0 else gauss2d_residual),
								array([0,-1,target_com[0],target_com[0],50,50,0]),
								args=(target[:,0],target[:,1],target[:,2]))
							params_frame.append(p_opt[0])
							target_zones_frame.append(target)
							maxhxys_frame.append(argmax([abs(gauss2dh(p_opt[0],i[0],i[1])) for i in target]))
							#---reverse curvature according to convention
							maxhs_frame.append(-1*gauss2dh(p_opt[0],target[maxhxys_frame[-1]][0],
								target[maxhxys_frame[-1]][1]))
						else:
							params_frame.append([])
							target_zones_frame.append([])
							maxhxys_frame.append([])
							maxhs_frame.append([])
					params.append(params_frame)
					maxhs.append(maxhs_frame)
					maxhxys.append(maxhxys_frame)
					target_zones.append(target_zones_frame)
				#---save the data to the MembraneData object
				result_data_collection = []
				for pnum in range(len(tl)):
					ps = tl[pnum]
					result_data = MembraneData('dimple2',label=sysname)
					for i in range(len(params)):
						#---store order is params, maxhs, maxhxys, target_zones
						result_data.add([params[i][pnum],maxhs[i][pnum],maxhxys[i][pnum],
							target_zones[i][pnum]],[range(len(mset.surf))[whichframes][i]])
					for i in analysis_descriptors[aname]:
						if i != 'testlist':
							result_data.addnote([i,(analysis_descriptors[aname])[i]])
					result_data.addnote(['testlist',testlist])	
					result_data.addnote(['this_test',pnum])	
					result_data.addnote(['cutoff',cutoff])
					result_data.addnote(['zfiltdir',zfiltdir])
					result_data.addnote(['decay_z0',decay_z0])
					result_data.addnote(['tl',tl])
					result_data_collection.append(result_data)
					del result_data
				for resdat in result_data_collection:
					msdats[anum].append(resdat)
				pickledump(msdats[anum],pklname,directory=pickles)
				del result_data_collection

#---plot summary
if 'plot' in routine:
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		msdat = unpickle_special('pkl.dimple2.'+filespec+'.pkl')
		anum = analysis_names.index(aname)
		msdats[anum] = msdat
	for param in params_plot:
		cutoff,zfiltdir = param
		#---lookup correct msdat items		
		msdats_subset = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			for ms in range(len(msdats[anum])):
				if msdats[anum][ms].getnote('cutoff') == cutoff and \
					msdats[anum][ms].getnote('zfiltdir') == zfiltdir:
					msdats_subset[anum].append(ms)
		#---plot
		fig = plt.figure(figsize=(10,len(analysis_names)*2))
		gs = gridspec.GridSpec(len(analysis_names),(3 if residplot == False else 4),wspace=0.0,hspace=0.0)
		axlist = []
		maxpeak = 0.
		#---limits of Hmax
		extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())])
			for mset in msets])
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			#---choose a test
			anum = analysis_names.index(aname)
			mset = msets[anum]
			dat = [msdats[anum][i] for i in msdats_subset[anum]]
			casual_names = [test.getnote('testlist')[test.getnote('this_test')] for test in dat]
			casual_names = [i[0] if type(i) == list else i for i in casual_names]
			color_order = [casual_names.index(i) for i in casual_names if i not in ['peak','valley','full']]
			#---plot Hmax		
			ax = plt.subplot(gs[anum,0:2])
			axlist.append(ax)
			if inset_raw_hmax or inset_extents: 
				axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,
					width="30%",height="30%",loc=1)
			ax.grid(True)
			if anum != len(analysis_names)-1: ax.set_xticklabels([])
			else: 
				ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
				plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
			ax.set_yticklabels([])
			ax.set_ylabel(label,fontsize=fsaxlabel)
			ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')
			for pnum in range(len(dat)):
				test = dat[pnum]
				if type(label) == list: label = label[0]
				if casual_names[pnum] == 'peak': color = colordict('blue')
				elif casual_names[pnum] == 'valley': color = colordict('red')
				elif casual_names[pnum] == 'full': color = colordict('black')
				else: color = colordict(color_order.index(pnum))
				
				#---find valid, fitted frames for calculating the mean
				fitted_inds = [type(test.data[i][1]) != list for i in range(len(test.data))]
####			#hmaxdat = array([i[1] for i in test.data if type(i[1]) != list])
				hmaxdat = array(test.data)[fitted_inds,1]
				#--->>>>>>>>?????????
				
				array(test.data)[:,1]				
				
				#---two-step filtering first to get valid fits and then to get this inside the curvature filter
				
				
				#cfilt_inds = list(where((1*array([type(test.data[i][1]) != list for i in range(len(test.data))])+1*array(abs(array(test.data)[:,1])>smallfilt)+1*array(abs(array(test.data)[:,1])<hifilt))==3)[0])

				cfilt_inds = [i for i in range(len(array(test.data))) if (type(array(test.data)[i][1]) != list and 10*abs(array(test.data)[i,1])>smallfilt and 10*abs(array(test.data)[i,1])<hifilt)]
				cfilt_inds_small = [i for i in range(len(array(test.data))) if (type(array(test.data)[i][1]) != list and 10*abs(array(test.data)[i,1])<smallfilt)]

				#cfilt_inds = list(where((1*array(fitted_inds)+1*array(abs(array(test.data)[:,1])>smallfilt)+1*array(abs(array(test.data)[:,1])<hifilt))==3)[0])
				
				hmaxdat = 10*array(test.data)[cfilt_inds,1]
				
				#---plot below small curvature limit for comparison
				hist,edges = numpy.histogram(array(test.data)[cfilt_inds_small,1],
					range=(-smallfilt,smallfilt),bins=1)
				ax.plot(1./2*(edges[1:]+edges[:-1]),hist,'o',c=color,
					lw=2)
				
				#---compute means
				#proper_subset = hmaxdat[array(abs(hmaxdat)>smallfilt)+array(abs(hmaxdat)<hifilt)]
				#valid_frames = sum((1*array(abs(hmaxdat)>smallfilt)+1*array(abs(hmaxdat)<hifilt))==2)
				valid_frames = len(hmaxdat)
				
				print 'pnum = '+str(pnum)+' mean Hmax = '+str(mean(hmaxdat))+\
					' valid = '+str('%1.2f'%valid_frames)
				#---filter by residuals
				if resid_filters != None and resid_filters != []:
					#unfilt_params = test.get(['type','params'])
					#valid_inds = [i for i in range(len(unfilt_params)) if type(unfilt_params[i]) != list]
					params = test.get(['type','params'])[cfilt_inds]
					maxhs = 10*test.get(['type','maxhs'])[cfilt_inds]
					target_zones = test.get(['type','target_zones'])[cfilt_inds]
					resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
						for j in range(len(target_zones))]
					means_by_resid_filt = []
					for i in range(len(resid_filters)):
						hist,edges = numpy.histogram(maxhs[array(resids)<resid_filters[i][0]],
							range=(-hifilt,hifilt),bins=hmax_nbins)				
						ax.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2,alpha=resid_filters[i][1])
						means_by_resid_filt.append(mean(maxhs[array(resids)<resid_filters[i][0]]))
				#---filtered data
				hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)>smallfilt],
					range=(-hifilt,hifilt),bins=hmax_nbins)
				if max(hist) > maxpeak: maxpeak = max(hist)
				ax.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,
					lw=2,label=(casual_names[pnum]+' ' if casual_names[pnum].find('monomer') == -1 else '')+\
					(str('%0.3f'%mean(hmaxdat)) if valid_frames > 0 else '0')+\
					(','+','.join([str('%0.3f'%i) 
						for i in means_by_resid_filt]) if (len(resid_filters) > 0 and valid_frames > 0) else '')+\
						r'$\,\mathrm{{nm}^{-1}}$('+\
						(str('%3.d'%(valid_frames)) if valid_frames > 0. else '0')+str(')'))
				#---inset with raw data, usually unnecessary if the fits are good enough
				#---Nb disabled this because most points were included			
				if inset_raw_hmax:
					hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)>hifilt],
						range=(-extrem_c0dat,extrem_c0dat),bins=hmax_nbins)
					axins.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2)
					hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)<hifilt],
						range=(-extrem_c0dat,extrem_c0dat),bins=hmax_nbins)
					axins.plot(1./2*(edges[1:]+edges[:-1]),hist,'.',c=color,lw=2)
				if inset_extents:
					extdat = sqrt(abs(array([test.data[i][0][4:6] for i in cfilt_inds \
						if type(test.data[i][1]) != list]).flatten()))
					hist,edges = numpy.histogram(extdat,
						range=(0,extent_range),bins=extent_range+1)
					axins.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2)
					'''
					if resid_filters != None and resid_filters != []:
						for i in range(len(resid_filters)):
							if len(extdat) > 0:
								hist,edges = numpy.histogram(extdat[array(resids)<resid_filters[i][0]],
									range=(0,extent_range),bins=extent_range+1)				
								axins.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2,
									alpha=resid_filters[i][1])
					'''
			ax.legend(loc='upper left',prop={'size':fsaxlegend_small})
			if inset_raw_hmax:
				axins.set_yticklabels([])
				axins.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=4))
				axins.set_xlabel('unfiltered')
			if inset_extents:
				axins.set_yticklabels([])
				axins.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=4))
				axins.set_xlabel('extent '+r'$\mathrm{\sigma_{a,b}\,(nm)}$')
			#---get protein
			if protein_pkl == None: mset_protein = msets[anum]
			#---retrieve protein points from a different dictionary item
			else:
				protein_pkl_name = protein_pkl
				for i in analysis_descriptors[protein_pkl_name]: 
					vars()[i] = (analysis_descriptors[protein_pkl_name])[i]
				filespec = specname_guess(sysname,trajsel)
				mset_protein = unpickle_special('pkl.structures.'+filespec+'.pkl')
				#---re-populate the variables for the original system now that you got the protein pickle
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			#---plot structure
			ax = plt.subplot(gs[anum,2])
			im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=mpl.cm.RdBu_r,
				inset=False,cmap_washout=1.0,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.])
			if analysis_names.index(aname) < len(analysis_names)-1:
				ax.set_yticklabels([])				
			#---height color scale
			axins2 = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins2,orientation="vertical")
			if inset_raw_hmax or inset_extents: plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins2.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',
				fontsize=fsaxlabel,rotation=270)
			#---plot protein hulls
			for pnum in range(len(dat)):
				test = dat[pnum]
				ps = test.getnote('tl')[test.getnote('this_test')]
				if type(ps) == slice:
					protpts = mean(mset_protein.protein,axis=0)[ps,0:2]
				else: 
					protpts = mean(mset_protein.protein,axis=0)[:,0:2]+ps
				casual_name = test.getnote('testlist')[test.getnote('this_test')]
				if casual_name != 'full' or len(casual_names) == 1:
					if casual_name == 'peak': color = colordict('blue')
					elif casual_name == 'valley': color = colordict('red')
					elif casual_name == 'full': color = colordict('black')
					else: color = colordict(color_order.index(pnum))
					#color='k'
					plothull(ax,protpts,mset=mset_protein,subdivide=None,c=color,alpha=1.)
				if pnum == 0 and residplot == True:
					axresid = plt.subplot(gs[anum,3])
					resid_scan = []
					for fr in range(len(target_zones)):
						H, xedges, yedges = histogram2d(
							linalg.norm(target_zones[fr][:,:2]-mean(protpts,axis=0),axis=1),
							[abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[fr]],
							range=((0,700),(0,10)),bins=(100,20))
						resid_scan.append(H)
					axresid.imshow(sum(resid_scan,axis=0).T,interpolation='nearest',origin='lower')
		for ax in axlist:
			ax.set_ylim((0,maxpeak))
		if zfiltdir == 0: filtstring = 'no'
		elif zfiltdir == -1: filtstring = 'down'
		elif zfiltdir == 1: filtstring = 'up'
		fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
		plt.savefig(pickles+'fig-dimple2-'+bigname+'.cut'+str(cutoff)+'.filt-'+\
			str(filtstring)+'.png',dpi=300,bbox_inches='tight')
		if show_plots: plt.show()
		plt.clf()
	
#---plot residual maps
if 'plot_residmaps' in routine:
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		msdat = unpickle_special('pkl.dimple2.'+filespec+'.pkl')
		anum = analysis_names.index(aname)
		msdats[anum] = msdat
	for param in params_plot:
		cutoff,zfiltdir = param
		#---lookup correct msdat items		
		msdats_subset = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			for ms in range(len(msdats[anum])):
				if msdats[anum][ms].getnote('cutoff') == cutoff and \
					msdats[anum][ms].getnote('zfiltdir') == zfiltdir:
					msdats_subset[anum].append(ms)
		#---plot
		fig = plt.figure(figsize=(10,len(analysis_names)*2))
		gs = gridspec.GridSpec(len(analysis_names),
			1+max([len([msdats[anum][i] for i in j]) for j in msdats_subset]),wspace=0.0,hspace=0.0)
		#---limits of Hmax
		extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())])
			for mset in msets])
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			#---choose a test
			anum = analysis_names.index(aname)
			mset = msets[anum]
			dat = [msdats[anum][i] for i in msdats_subset[anum]]
			casual_names = [test.getnote('testlist')[test.getnote('this_test')] for test in dat]
			casual_names = [i[0] if type(i) == list else i for i in casual_names]
			color_order = [casual_names.index(i) for i in casual_names if i not in ['peak','valley','full']]
			for pnum in range(len(dat)):
				test = dat[pnum]
				if type(label) == list: label = label[0]
				if casual_names[pnum] == 'peak': color = colordict('blue')
				elif casual_names[pnum] == 'valley': color = colordict('red')
				elif casual_names[pnum] == 'full': color = colordict('black')
				else: color = colordict(color_order.index(pnum))
				casual_name = test.getnote('testlist')[test.getnote('this_test')]
				#---plot the residual maps for individual proteins
				casual_name = test.getnote('testlist')[test.getnote('this_test')]
				unfilt_params = test.get(['type','params'])
				valid_inds = [i for i in range(len(unfilt_params)) if type(unfilt_params[i]) != list]
				params = test.get(['type','params'])[valid_inds]
				maxhs = test.get(['type','maxhs'])[valid_inds]
				target_zones = test.get(['type','target_zones'])[valid_inds]
				ax = plt.subplot(gs[anum,1+pnum])
				print 'status: computing residual map, '+str(aname)+', '+str(pnum)
				vecs = mean(mset.vecs,axis=0)
				m,n = mset.griddims
				getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
					for i in linspace(0,vecs[0]/mset.lenscale,m)])
				sum_resids = [[[] for i in range(n+10)] for j in range(m+10)]
				for fr in range(len(target_zones)):
					gridpts = [[int(i[d]/mset.vec(fr)[d]*(mset.griddims[d]-1)) for d in range(2)] 
						for i in target_zones[fr][:,:2]]
					resids = [i[2]-gauss2d(params[fr],i[0],i[1]) for i in target_zones[fr]]
					for k in range(len(gridpts)):
						sum_resids[gridpts[k][0]][gridpts[k][1]].append(resids[k])
				sum_resids = [[sum(sum_resids[j][i]) for i in range(n)] for j in range(m)]
				extrem = max([abs(array(sum_resids).min()),array(sum_resids).max()])
				ax.imshow(array(array(sum_resids)[:m,:n]).T,interpolation='nearest',origin='lower',
					cmap=mpl.cm.RdBu_r,vmin=-extrem,vmax=extrem)
				ax.set_yticklabels([])
				ax.set_xticklabels([])
			#---get protein
			if protein_pkl == None: mset_protein = msets[anum]
			#---retrieve protein points from a different dictionary item
			else:
				protein_pkl_name = protein_pkl
				for i in analysis_descriptors[protein_pkl_name]: 
					vars()[i] = (analysis_descriptors[protein_pkl_name])[i]
				filespec = specname_guess(sysname,trajsel)
				mset_protein = unpickle_special('pkl.structures.'+filespec+'.pkl')
				#---re-populate the variables for the original system now that you got the protein pickle
				for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			#---plot structure
			ax = plt.subplot(gs[anum,0])
			im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,lognorm=False,cmap=mpl.cm.RdBu_r,
				inset=False,cmap_washout=1.0,ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.])
			ax.set_yticklabels([])
			ax.set_xticklabels([])
			ax.set_xlabel('')
			ax.set_ylabel(label)
			#---plot protein hulls
			for pnum in range(len(dat)):
				test = dat[pnum]
				ps = test.getnote('tl')[test.getnote('this_test')]
				if type(ps) == slice:
					protpts = mean(mset_protein.protein,axis=0)[ps,0:2]
				else: 
					protpts = mean(mset_protein.protein,axis=0)[:,0:2]+ps
				casual_name = test.getnote('testlist')[test.getnote('this_test')]
				if casual_name != 'full' or len(casual_names) == 1:
					if casual_name == 'peak': color = colordict('blue')
					elif casual_name == 'valley': color = colordict('red')
					elif casual_name == 'full': color = colordict('black')
					else: color = colordict(color_order.index(pnum))
					plothull(ax,protpts,mset=mset_protein,subdivide=None,c=color,alpha=1.)
				if pnum == 0 and residplot == True:
					axresid = plt.subplot(gs[anum,3])
					resid_scan = []
					for fr in range(len(target_zones)):
						H, xedges, yedges = histogram2d(
							linalg.norm(target_zones[fr][:,:2]-mean(protpts,axis=0),axis=1),
							[abs(gauss2d(params[0],i[0],i[1])-i[2])**2 for i in target_zones[fr]],
							range=((0,700),(0,10)),bins=(100,20))
						resid_scan.append(H)
					axresid.imshow(sum(resid_scan,axis=0).T,interpolation='nearest',origin='lower')
		if zfiltdir == 0: filtstring = 'no'
		elif zfiltdir == -1: filtstring = 'down'
		elif zfiltdir == 1: filtstring = 'up'
		fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
		plt.savefig(pickles+'fig-dimple2-'+bigname+'-residual_map.cut'+str(cutoff)+'.filt-'+\
			str(filtstring)+'.png',dpi=300,bbox_inches='tight')
		if show_plots: plt.show()
		plt.clf()


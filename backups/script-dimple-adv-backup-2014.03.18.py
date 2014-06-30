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
		'trajsel':'s6-sim-lonestar/md.part0002.400000-140000-200.xtc',
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
		'label':r'$\mathrm{control\,(v2)}$',
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
routine = ['compute','plot'][1:]
bigname = 'v614-v612-v550-combine'

#---parameter set for sweeps
params = [[c,d] for c in [50,80,100,120,150,200,250] for d in [0,1,-1]]

#---parameters
cutoff = 320
decay_z0 = False

#---plot settings
hifilt = 0.05
smallfilt = 0.001
extrem_c0_limit = 10.

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
if 'msets' not in globals() and 'plot' not in routine:
	msets = [[] for i in range(len(analysis_names))]
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		msdat = unpickle_special('pkl.dimple2.'+filespec+'.cut'+str(cutoff)+'.pkl')
		anum = analysis_names.index(aname)
		mset = unpickle_special('pkl.structures.'+filespec+'.pkl')
		msets[anum] = mset
		#msdats[anum] = msdat

#---fits
if 'compute' in routine:
	for cutoff in [50,80,100,120,150,200,250]:
		for aname in analysis_names:
			result_data_collection = []
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---only run the computation if the pkl was missing during the load step
			if msdats[anum] == [] or msdats[anum] == None:
				filespec = specname_guess(sysname,trajsel)
				pklname = 'pkl.dimple2.'+filespec+'.cut'+str(cutoff)+'.pkl'
				print 'status: running computation to generate '+str(pklname)
				#---protein accounting
				mset = msets[anum]
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
				#---prepare a list of tests
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
				#---perform fits over protein slices
				params = []
				maxhs = []
				maxhxys = []
				maxhxys = []
				target_zones = []
				for fr in range(len(mset.surf))[whichframes]:
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
						#---find the center of the target as initial guess for the fitter
						target_com = [mean(target[:,0]),mean(target[:,1])]
						#---perform the fit
						#---recall z0,c0,x0,y0,sx,sy,th = params
						p_opt = leastsq((gauss2d_residual_z0 if decay_z0 else gauss2d_residual),
							array([0,-1,target_com[0],target_com[0],50,50,0]),
							args=(target[:,0],target[:,1],target[:,2]))
						params_frame.append(p_opt[0])
						target_zones_frame.append(target)
						maxhxys_frame.append(argmax([abs(gauss2dh(p_opt[0],i[0],i[1])) for i in target]))
						#---reverse curvature according to convention
						maxhs_frame.append(-1*gauss2dh(p_opt[0],target[maxhxys_frame[-1]][0],
							target[maxhxys_frame[-1]][1]))
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
						result_data.add([params[i][pnum],maxhs[i][pnum],maxhxys[i][pnum],target_zones[i][pnum]],
							[range(len(mset.surf))[whichframes][i]])
					for i in analysis_descriptors[aname]:
						if i != 'testlist':
							result_data.addnote([i,(analysis_descriptors[aname])[i]])
					result_data.addnote(['testlist',testlist])	
					result_data.addnote(['this_test',pnum])	
					result_data.addnote(['cutoff',cutoff])
					result_data.addnote(['decay_z0',cutoff])
					result_data.addnote(['tl',tl])
					result_data_collection.append(result_data)
					del result_data
				pickledump(result_data_collection,pklname,directory=pickles)
				#del result_data_collection

# latpos = array([i[3][i[2]][:2] for i in test.data])
# plt.scatter(latpos[:,0],latpos[:,1]);plt.show()
 
#---plot summary
if 'plot' in routine:
	for cutoff in [100]:
		#---load for plots, specifically, if they are not already loaded
		if 'msets' not in globals():
			msets = [[] for i in range(len(analysis_names))]
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			filespec = specname_guess(sysname,trajsel)
			msdat = unpickle_special('pkl.dimple2.'+filespec+'.cut'+str(cutoff)+'.pkl')
			anum = analysis_names.index(aname)
			if msets[anum] == [] or msets[anum] == None:
				mset = unpickle_special('pkl.structures.'+filespec+'.pkl')
				msets[anum] = mset	
			msdats[anum] = msdat
		#---plot
		fig = plt.figure(figsize=(10,len(analysis_names)*2))
		gs = gridspec.GridSpec(len(analysis_names),3,wspace=0.0,hspace=0.0)
		#---limits of Hmax
		extrem_c0dat = max([max([[abs(j.max()),abs(j.min())] for j in [array([[i[0][1] 
			for i in test.data] for test in dat])]][0]) for dat in msdats])
		if extrem_c0dat > extrem_c0_limit: extrem_c0dat = extrem_c0_limit
		extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())]) for mset in msets])
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			#---choose a test
			anum = analysis_names.index(aname)
			mset = msets[anum]
			dat = msdats[anum]
			casual_names = [test.getnote('testlist')[test.getnote('this_test')] for test in dat]
			casual_names = [i[0] if type(i) == list else i for i in casual_names]
			color_order = [casual_names.index(i) for i in casual_names if i not in ['peak','valley','full']]
			#---plot Hmax		
			ax = plt.subplot(gs[anum,0:2])
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="30%",height="30%",loc=2)
			ax.grid(True)
			if anum != len(analysis_names)-1: ax.set_xticklabels([])
			else: 
				ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
				plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
			ax.set_yticklabels([])
			ax.set_ylabel(label,fontsize=fsaxlabel)
			for pnum in range(len(dat)):
				test = dat[pnum]
				if type(label) == list: label = label[0]
				if casual_names[pnum] == 'peak': color = colordict('blue')
				elif casual_names[pnum] == 'valley': color = colordict('red')
				elif casual_names[pnum] == 'full': color = colordict('black')
				else: color = colordict(color_order.index(pnum))
				#---plot below small curvature limit for comparison
				hmaxdat = array([i[1] for i in test.data])
				hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)<smallfilt],range=(-smallfilt,smallfilt),bins=21)
				ax.plot(1./2*(edges[1:]+edges[:-1]),hist,':',c=color,
					lw=2)
				#---filtered data
				hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)>smallfilt],range=(-hifilt,hifilt),bins=21)
				ax.plot(1./2*(edges[1:]+edges[:-1]),hist,c=color,
					lw=2,label=casual_names[pnum])
				#---inset with raw data
				c0dat = array([i[0][1] for i in test.data])
				hist,edges = numpy.histogram(c0dat,range=(-extrem_c0dat,extrem_c0dat),bins=21)
				axins.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2)
			ax.legend()
			axins.set_yticklabels([])
			axins.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=4))
			axins.set_xlabel('unfiltered')
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
			#---height color scale
			axins2 = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins2,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
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
					plothull(ax,protpts,mset=mset_protein,subdivide=None,c=color,alpha=1.)
		fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
		plt.savefig(pickles+'fig-dimple2-'+bigname+'.cut'+str(cutoff)+'.png',dpi=500,bbox_inches='tight')
		#plt.show()
		plt.clf()
		del msdats
		del msets
	

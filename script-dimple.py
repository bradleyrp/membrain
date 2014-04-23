#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

from scipy.optimize import leastsq

'''
Dimple-fitting code
This is the latest version created 2014.04.20
Here I introduce "dimple3" pickle objects and try to make the bookkeeping easier.

Code plan. 
	Make the membrane data object
	Make a function which designs the note objects
	Make a function which checks for pkls of previously created parameter sets
	Check the exponent on the fitting
	Check the RMSD
	Make a lattice code for testing the control
'''

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---additions to the library of available simulations
analysis_descriptors_extra = {
	'v614-120000-220000-200': {},
	'v612-75000-175000-200': {},
	'v550-400000-500000-160': {},
	'v614-40000-140000-200': {},
	'v612-10000-80000-200': {},
	'v550-300000-400000-200': {},
	'v616-210000-310000-200': {},
	}
		
#---coallate two dictionaries
execfile('header-cgmd.py')
master_dict = dict()
for key in analysis_descriptors.keys():
	if key in analysis_descriptors_extra.keys():
		master_dict.update({key:dict(analysis_descriptors[key],
			**analysis_descriptors_extra[key])})
	else: master_dict.update({key:analysis_descriptors[key]})
analysis_descriptors = master_dict

def dimple_compile_notes(aname,expt):
	'''Compile a "notes" object with all parameters for a test+system combination.'''
	notes_list = []
	for i in analysis_descriptors[aname]:
		notes_list.append([i,(analysis_descriptors[aname])[i]])
	for key in expt.keys():
		notes_list.append([key,expt[key]])
	notes_list.append(['callsign',aname])
	return notes_list

#---analysis plan
analysis_names = [
	'v616-210000-310000-200',
	'v614-120000-220000-200',
	'v612-75000-175000-200',
	'v550-400000-500000-160',
	'v614-40000-140000-200',
	'v612-10000-80000-200',
	'v550-300000-400000-200',
	][:4]
routine = [
	'compute',
	'plot',
	][1:]
bigname = 'v616-v614-v612-v550-ver2'

#---EXPERIMENT PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---global settings
gridspacing = 1.0

'''
Experiment/analysis parameters list
	1. decay_z0 : 
		zero = dimple decays to zero (average height) at infinity
		min = dimple decays to the average bilayer minimum height at infinity
	2. cutoff : cutoff about protein points which defines a neighborhood (in Angstroms)
	3. z_filter : 
		inclusive = 
		up =
		down =
	4. fit_dynamics_type : 
		dynamic = 
		mean = 
	5. framecounts : if None use all frames, if slice or integer, use that to get a subset of frames
	6. geography : 
		proteins = use the proteins themselves as geographic sources
		lattice = sample random points and use a fixed single point as the center of the neighborhood
'''

#---experiment parameters
params_expt = {
	'decay_z0' : 'zero',
	'cutoff' : 100,
	'fit_dynamics_type' : 'dynamic',
	'z_filter' : 'inclusive',
	'framecounts' : 500, 
	'geography' : '',
	}
	
#---sweeping parameters
params_sweeps = {
	'cutoff':[100],
	'geography':[
		'proteins',
		['control','v614-120000-220000-200'],
		['control','v612-75000-175000-200'],
		'max','min',
		'max_dynamic','min_dynamic',
		['lattice_grid',6],
		][:],
	}	
	
#---PLOT PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---general parameters
show_plots = True

#---general settings
params_plot_settings = {
	'hifilt' : 0.06,
	'smallfilt' : 0.001,
	'hist_step' : 0.005,
	'extent_range' : 32,
	}

#---master list of parameters
params_plot_master,params_plot_master_names = [],[]

#---plot parameter set "lattice6"
params_plot_master_names.append('lattice6')
params_plot_master.append([
	[dimple_compile_notes(aname,
	dict({
		'decay_z0' : 'zero',
		'cutoff' : 100,
		'fit_dynamics_type' : 'dynamic',
		'z_filter' : 'inclusive',
		'framecounts' : 500, 
		'geography' : geog,
		},))
		for geog in [['lattice_grid',6]]]
	for aname in analysis_names])

#---plot parameter set "proteins_controls"
params_plot_master_names.append('proteins_controls')
params_plot_master.append(
	[[dimple_compile_notes(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			},))
		for geog in ['proteins']] 
			for aname in analysis_names 
			if (analysis_descriptors[aname])['nprots'] > 0]+\
	[[dimple_compile_notes(aname,
	dict({
		'decay_z0' : 'zero',
		'cutoff' : 100,
		'fit_dynamics_type' : 'dynamic',
		'z_filter' : 'inclusive',
		'framecounts' : 500, 
		'geography' : geog,
		},))
	for geog in [['control','v614-120000-220000-200'],['control','v612-75000-175000-200']]
		for aname in analysis_names
		if (analysis_descriptors[aname])['nprots'] == 0]])

#---plot parameter set "extrema"
params_plot_master_names.append('extrema')
params_plot_master.append(
	[[dimple_compile_notes(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			},))
		for geog in ['max_dynamic','min_dynamic','max','min']] 
			for aname in analysis_names])

#---select a plot scheme
params_plot = params_plot_master[0]
plot_scheme_name = params_plot_master_names[params_plot_master.index(params_plot)]

#---NOTES
#-------------------------------------------------------------------------------------------------------------

'''
COMPUTE WORKFLOW
	1. Set experiment or analysis parameters and choose systems.
	2. If sweeping, make a set of dictionaries with the parameters.
	3. Over each sub-experiment, get the systems and assemble the notes list.
	4. Use the notes list to see if you have already done that experiment.
	5. If not, run the experiment and store the notes list.
	6. Save the pickle either after each experiment, or at the end.
'''
	
#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---NOTE
#---The function gauss2dh returns a value equal to 2H, which is corrected with the following factor.
curvfac = float32(0.5)

#---load expressions for mean and Gaussian curvature
execfile('script-curvature-calculus.py')

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)
def gauss2d_residual(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual.'''
	return ((gauss2d(params,x,y)-z)**2)
	
def gauss2d_z0(params,x,y,z0=0):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, decay height fixed.'''
	c0,x0,y0,sx,sy,th = params[1:]
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)
def gauss2d_residual_z0(params,x,y,z,z0=0):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual, decay height fixed'''
	return ((gauss2d_z0(params,x,y,z0=z0)-z))

def gauss2d_z0_global(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, decay height fixed.'''
	c0,x0,y0,sx,sy,th = params[1:]
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)
def gauss2d_residual_z0_global(params,x,y,z):
	'''Two-dimensional Gaussian height function with fluid axis of rotation, residual, decay height fixed'''
	return ((gauss2d_z0(params,x,y,z0=z0)-z))

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---tag for referring to gridspacing/rounder in structure objects
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'

def dimple_save_data(aname):
	'''Saves the global dimple3 result data object for a particular system name (name msdats[anum]).'''
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	anum = analysis_names.index(aname)
	filespec = specname_guess(sysname,trajsel)
	pklname = 'pkl.dimple3.'+spacetag+filespec+'.pkl'
	pickledump(msdats[anum],pklname,directory=pickles)

def listlook(data,name):
	'''Treat a list like a dictionary and look up a value.'''
	for i in data:
		if i[0] == name:
			return i[1]

def dimple_test_lookup(notes,aname):
	'''Find the indices corresponding to a particular test.'''
	anum = analysis_names.index(aname)
	inds = [j for j in range(len(msdats[anum])) if all([i in msdats[anum][j].notes for i in notes]) == True]
	if inds != []: return inds
	else: return False

def dimple_generate_neighborhood(notes):
	'''Generate neighborhood points.'''
	request_geography = listlook(notes,'geography')
	#---select frames
	framecounts = listlook(notes,'framecounts')
	if framecounts != None and type(framecounts) != slice: frame_select = range(framecounts)
	elif type(framecounts) == slice: frame_select = range(len(mset_protein.protein))[framecounts]
	else: frame_select = range(len(mset_protein.protein))
	#---select the surface
	aname = listlook(notes,'callsign')
	anum = analysis_names.index(aname)
	mset = msets[anum]
	#---geography is based on the protein positions
	if request_geography == 'proteins' or request_geography[0] == 'control':
		#---get protein points
		if request_geography[0] == 'control':
			print 'status: test = control, system = '+str(request_geography[1])
			callsign_protein = listlook(notes,'geography')[1]
			sysname,trajsel = [(analysis_descriptors[callsign_protein])[key] 
				for key in ['sysname','trajsel']]
			filespec = specname_guess(sysname,trajsel)
			mset_protein = unpickle(pickles+'pkl.structures.'+spacetag+filespec+'.pkl')
			nprots = (analysis_descriptors[callsign_protein])['nprots']
		else: 
			mset_protein = mset
			nprots = listlook(notes,'nprots')
		#---send protein neighborhoods based on monomer definitions
		if nprots == 1 or request_geography[0] == 'control':
			print 'status: test = proteins, nprots = 1'
			nameslist = ['monomer']
			print nameslist
			return [[mset_protein.protein[fr][:,0:2] for fr in frame_select]],nameslist
		elif nprots > 1:
			print 'status: test = proteins, nprots = '+str(nprots)
			nresprot = shape(mset_protein.protein)[1]/nprots
			prot_slices = [slice(i*nresprot,(i+1)*nresprot) for i in range(nprots)]
			collected_pts = [[mset_protein.protein[fr][:,0:2] for fr in frame_select]]
			for sl in prot_slices:
				collected_pts.append([mset_protein.protein[fr][sl,0:2] for fr in frame_select])
			nameslist = ['oligomer']+['monomer '+str(i+1) for i in range(len(prot_slices))]
			return collected_pts,nameslist
	#---use the average midplane maximum height as a fixed point
	elif request_geography == 'max':
		maxpos_disc = unravel_index(mean(mset.surf,axis=0).argmax(),mset.griddims)
		maxpos = [float(maxpos_disc[i])/mset.griddims[i]*mean(mset.vecs,axis=0)[i] for i in range(2)]
		return [[[maxpos] for fr in frame_select]],['peak']
	#---use the average midplane minimum height as a fixed point
	elif request_geography == 'min':
		minpos_disc = unravel_index(mean(mset.surf,axis=0).argmin(),mset.griddims)
		minpos = [float(minpos_disc[i])/mset.griddims[i]*mean(mset.vecs,axis=0)[i] for i in range(2)]
		return [[[minpos] for fr in frame_select]],['valley']
	#---use the instantaneous midplane maximum height as a fixed point
	elif request_geography == 'max_dynamic':
		maxpos_disc = [unravel_index(mset.surf[fr].argmax(),mset.griddims) for fr in frame_select]
		return [[[[float(maxpos_disc[fr][i])/mset.griddims[i]*mset.vec(fr)[i] 
			for i in range(2)]] for fr in frame_select]],['peak (dynamic)']
	#---use the instantaneous midplane maximum height as a fixed point
	elif request_geography == 'min_dynamic':
		minpos_disc = [unravel_index(mset.surf[fr].argmin(),mset.griddims) for fr in frame_select]
		return [[[[float(minpos_disc[fr][i])/mset.griddims[i]*mset.vec(fr)[i] 
			for i in range(2)]] for fr in frame_select]],['valley (dynamic)']
	#---choose static points according to a lattice grid with a specific number of points per dimension
	elif request_geography[0] == 'lattice_grid':
		latpts = [[[[i,j]] for j in linspace(0,mset.vec(fr)[1],request_geography[1]*2+1)[1::2] 
			for i in linspace(0,mset.vec(fr)[0],request_geography[1]*2+1)[1::2]] for fr in frame_select]
		return array(latpts).transpose((1,0,2,3)),['position '+str(i+1) for i in range(len(latpts[0]))]
	else: print 'error: cannot understand neighborhood request'	
		
def dimple_colordict(ask,listing=None,scheme=None,nprots=None):
	'''Color lookup for plotting protein hulls.'''
	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in [0,1,2,3,4,6,7,8]]
	clrsdict_specific = {'red':clrs[0],'blue':clrs[1],'black':'k'}
	clrsdict_leftovers = list(array(clrs)[range(2,len(clrs))])
	if scheme != None and re.match('lattice.+',scheme):
		clrsdict_leftovers = [brewer2mpl.get_map('Set3','qualitative',12).mpl_colors[i] for i in range(12)]
	reserved_colors = ['max','min,','max (dynamic)','min (dynamic)','oligomer']
	if type(ask) == str:
		if ask == 'peak' or ask == 'peak (dynamic)': return clrsdict_specific['red']
		elif ask == 'valley' or ask == 'valley (dynamic)': return clrsdict_specific['blue']
		elif ask == 'oligomer' or nprots == 1: return clrsdict_specific['black']
		elif listing != None:
			nonreserved_names = [i for i in listing if i not in reserved_colors]
			if ask in nonreserved_names:
				return clrsdict_leftovers[nonreserved_names.index(ask)%len(clrsdict_leftovers)]
		else: return 'y'
	elif type(ask) == int:
		return clrsdict_leftovers[ask%len(clrsdict_leftovers)]		
	else: return 'y'
	
def dimple_colordict_systems(nname):
	'''This function handles colors and names for tests of the individual systems.'''
	bold_blue = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[1]
	bold_green = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[2]
	bold_purple = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[3]
	bold_light_purple = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[7]
	bold_gray = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[7]
	bold_black = 'k'
	colornamedict = {
		'8xENTH(close)' : [bold_blue,r'$\mathrm{8 \times ENTH}$'],
		'4xENTH' : [bold_green,r'$\mathrm{4 \times ENTH}$'],
		'1xENTH' : [bold_light_purple,r'$\mathrm{1 \times ENTH}$'],
		'control' : [bold_black,r'$\mathrm{control}$'],
		}
	return colornamedict[nname]
	
#---MAIN, COMPUTE
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		filespec = specname_guess(sysname,trajsel)
		anum = analysis_names.index(aname)
		mset = unpickle(pickles+'pkl.structures.'+spacetag+filespec+'.pkl')
		msets[anum] = mset

#---compute fits
if 'compute' in routine:

	#---load the current set of msdats which may already have some of the requested tests in the sweep
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		#---load the pickle to see if it exists	
		filespec = specname_guess(sysname,trajsel)
		pklname = 'pkl.dimple3.'+spacetag+filespec+'.pkl'
		msdat = unpickle(pickles+pklname)
		if msdat == None: msdats[anum] = []
		else: msdats[anum] = msdat

	#---prepare a list of experiments according to any requested sweeps
	list_params_expts = []
	keylist = params_sweeps.keys()
	sweeps = [params_sweeps[key] for key in keylist]
	param_combos = list(itertools.product(*sweeps))
	for new_expt in param_combos:
		expt_copy = params_expt.copy()
		for key in keylist:
			expt_copy[key] = new_expt[keylist.index(key)]
		list_params_expts.append(expt_copy)

	#---loop over experiments and systems and combine valid experiments into a list of notes objects
	list_expts = []
	for expt in list_params_expts:
		for aname in analysis_names:
			notes = dimple_compile_notes(aname,expt)
			#---remove contradictory protein/control tests and see if test was already completed
			if not (listlook(notes,'nprots') == 0 and listlook(notes,'geography') == 'proteins') and \
				not (listlook(notes,'nprots') > 0 and listlook(notes,'geography')[0] == 'control'):				
					list_expts.append(notes)
			
	#---master loop over all dimple fits
	for notes in list_expts:
		aname = listlook(notes,'callsign')
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		mset = msets[anum]
		
		#---define the neighborhood
		print '\n',
		nborhoods,nbornames = dimple_generate_neighborhood(notes)
		print 'status: neighborhood shape = '+str(shape(nborhoods))
		print 'status: neighborhood names = '+str(nbornames)
		
		#---select frames
		framecounts = listlook(notes,'framecounts')
		if framecounts != None and type(framecounts) != slice: frame_select = range(framecounts)
		elif type(framecounts) == slice: frame_select = range(len(mset_protein.protein))[framecounts]
		else: frame_select = range(len(mset_protein.protein))

		#---loop over neighborhoods within each test type
		for nborhood_ind in range(len(nborhoods)):
			nborhood = nborhoods[nborhood_ind]
			
			#---check to see if this calculation was already completed
			notes_nbor = notes+[['neighborhood_name',nbornames[nborhood_ind]]]
			if (False == dimple_test_lookup(notes_nbor,aname)):
			
				#---prepare the result file
				result_data = MembraneData('dimple3',label=sysname)
				result_data.notes = notes_nbor
					
				#---unpack the notes
				for i in notes_nbor: vars()[i[0]] = i[1]
		
				#---dimple fitting code
				print 'status: running test'
				params = []
				maxhs = []
				maxhxys = []
				maxhxys = []
				target_zones = []
				for fr in frame_select:
					status('neighborhood = '+str(nborhood_ind+1)+'/'+str(len(nborhoods))+' frame = '+str(fr))
					#---replicate surface points under periodic boundary conditions
					surfpts = mset.wrappbc(mset.unzipgrid(mset.surf[fr],vecs=mset.vec(0)),
						vecs=mset.vec(fr),mode='nine')
					#---find minimum distances
					cd = scipy.spatial.distance.cdist(nborhood[fr],surfpts[:,:2])
					tmp = array([np.min(cd,axis=0),surfpts[:,2]]).T
					selected = where(tmp[:,0]<cutoff)[0]
					target = surfpts[selected]
					if z_filter == 'up':
						target = target[target[:,2]>0]
					elif z_filter == 'down':
						target = target[target[:,2]<0]
					elif z_filter == 'inclusive': 1				
					else: print 'error: cannot understand filter'
					#---find the center of the target as initial guess for the fitter
					target_com = [mean(target[:,0]),mean(target[:,1])]
					#---perform the fit
					#---recall z0,c0,x0,y0,sx,sy,th = params
					#---check for a sufficient number of points in the neighborhood
					if len(target) >= 7:
						if decay_z0 == 'zero': residfunc = gauss2d_residual_z0
						elif decay_z0 == 'min':
							residfunc = gauss2d_residual_z0_global
							z0 = mean(mset.surf,axis=0).min()
							print 'status: set decay to average bilayer height minimum'
						else: print 'error: cannot understand decay_z0 setting'
						p_opt = leastsq(residfunc,array([0,0,target_com[0],target_com[0],50,50,0]),
							args=(target[:,0],target[:,1],target[:,2]))
						params.append(p_opt[0])
						target_zones.append(target)
						maxhxys.append(argmax([abs(gauss2dh(p_opt[0],i[0],i[1])) for i in target]))
						#---reverse curvature according to convention
						maxhs.append(-1*gauss2dh(p_opt[0],target[maxhxys[-1]][0],
							target[maxhxys[-1]][1]))
					else:
						params.append([])
						target_zones.append([])
						maxhxys.append([])
						maxhs.append([])

				#---populate the result object with the dimple fit results
				for i in range(len(params)):
					#---store order is params, maxhs, maxhxys, target_zones
					result_data.add([params[i],maxhs[i],maxhxys[i],target_zones[i]],[frame_select[i]])
				msdats[anum].append(result_data)
				del result_data
				print '\n',
				dimple_save_data(aname)

#---MAIN, PLOT
#-------------------------------------------------------------------------------------------------------------

#---plot summary
if 'plot' in routine:
	
	if 'msdats' not in globals():
		#---load the current set of msdats
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---load the pickle to see if it exists	
			filespec = specname_guess(sysname,trajsel)
			pklname = 'pkl.dimple3.'+spacetag+filespec+'.pkl'
			msdat = unpickle(pickles+pklname)
			if msdat == None: msdats[anum] = []
			else: msdats[anum] = msdat
			
	#---prepare plot panels	
	fig = plt.figure(figsize=(10,len(params_plot)*2))
	gs = gridspec.GridSpec(len(params_plot),3,wspace=0.0,hspace=0.0)
	gs.update(left=0.0,right=0.7)
	gs2 = gridspec.GridSpec(len(params_plot),1,wspace=0.0,hspace=0.0)
	gs2.update(left=0.75,right=1.0)
	axlist,axlist_struct,axlist_extent = [],[],[]
	fp = open(pickles+'calc-dimple3-'+bigname+\
		'.cut'+str(listlook(params_plot[0][0],'cutoff'))+spacetag+\
		plot_scheme_name+'.txt','w')

	#---global plot specifications
	extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())])
		for mset in msets])
	maxpeak,maxpeak_extent = 0.,0.
	phists = []
	lattice_mean_collect = [[] for j in range(len(params_plot))]

	#---loop over rows
	for pnum in range(len(params_plot)):
		ax = plt.subplot(gs[pnum,0:2])
		axlist.append(ax)
		ax_struct = plt.subplot(gs2[pnum])
		axlist_struct.append(ax_struct)
		ax_extent = plt.subplot(gs[pnum,2])
		axlist_extent.append(ax_extent)
		nnames_control = []
		
		#---loop over plots within a panel
		chists = []
		for gnum in range(len(params_plot[pnum])):
			expt = params_plot[pnum][gnum]
			aname = listlook(expt,'callsign')
			anum = analysis_names.index(aname)
			mset = msets[anum]
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			nborhoods,nbornames = dimple_generate_neighborhood(expt)
			#---write header for the calculations text output
			fp.write('\nsystem = '+str(aname)+'\t'+(analysis_descriptors[aname])['label_text']+'\n\n')
			table_top = 'name        '+'<H_max> (nm^-1)   '+'frames  '+'RMSD(A) '+\
				'sigma_a,b (nm)  '+'sigma_a,b unfiltered (nm)'
			fp.write(table_top+'\n')
			fp.write('-'.join(['' for i in range(len(table_top)+1)])+'\n')
	
			#---get indices for plotting
			dat_inds = dimple_test_lookup(expt,aname)
			#---retreive names for color selection
			nnames = [listlook(msdats[anum][cnum].notes,'neighborhood_name') for cnum in dat_inds]
			#---get general plot settings
			for i in params_plot_settings: vars()[i] = params_plot_settings[i]
			hmax_nbins = int(2*hifilt/hist_step)+1
			
			#---loop over curves within one plot
			plot_order = range(len(dat_inds))
			if 'oligomer' in nbornames: plot_order = range(len(dat_inds))[1:]+[0]
			for cnum in plot_order:
				ax = axlist[pnum]
				test = msdats[anum][dat_inds[cnum]]
				#---set legend labels and colors for control
				geog = listlook(test.notes,'geography')
				if geog[0] == 'control':
					label = dimple_colordict_systems((analysis_descriptors[geog[1]])['label_text'])[1]
					nnames_control.append(label)
					#---previously used an ad hoc color set
					if 0: color = dimple_colordict(label,listing=nnames_control)
					#---use the color dedicated to the protein system the control used as hypothesis
					color = dimple_colordict_systems((analysis_descriptors[geog[1]])['label_text'])[0]
				else: 
					label = listlook(test.notes,'neighborhood_name')
					if listlook(test.notes,'nprots') == 1 and plot_scheme_name != 'extrema' and \
						re.match('lattice.+',plot_scheme_name) == None:
						nprots = 1
					else: nprots = None
					nprots = (listlook(test.notes,'nprots') if \
							re.match('lattice.+',plot_scheme_name) == None else None)
					color = dimple_colordict(nnames[cnum],listing=nnames,
						scheme=plot_scheme_name,nprots=nprots)
				if label in ['peak (dynamic)','valley (dynamic)']: alpha = 0.5
				else: alpha = 1.
				status('panel = '+str(pnum+1)+'/'+str(len(params_plot))+
					' curve = '+str(cnum+1)+'/'+str(len(dat_inds)))
			
				#---filtering
				fitted_inds = [type(test.data[i][1]) != list for i in range(len(test.data))]
				hmaxdat = curvfac*array(test.data)[fitted_inds,1]
				#---two-step filtering first to get valid fits and then to apply the curvature filter
				cfilt_inds = [i for i in range(len(array(test.data))) 
					if (type(array(test.data)[i][1]) != list 
					and curvfac*10*abs(array(test.data)[i,1])>smallfilt 
					and curvfac*10*abs(array(test.data)[i,1])<hifilt)]
				hmaxdat = curvfac*10*array(test.data)[cfilt_inds,1]
				#---report residuals
				params = test.get(['type','params'])[cfilt_inds]
				maxhs = curvfac*10*test.get(['type','maxhs'])[cfilt_inds]
				target_zones = test.get(['type','target_zones'])[cfilt_inds]
				resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 for i in target_zones[j]])) 
					for j in range(len(target_zones))]
				#---filtered data
				hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)>smallfilt],
					range=(-hifilt-hist_step/2,hifilt+hist_step/2),bins=hmax_nbins)
				if max(hist) > maxpeak: maxpeak = max(hist)
				ax.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2,label=label,alpha=alpha)
				hmax_mids = 1./2*(edges[1:]+edges[:-1])
				if pnum == len(params_plot)-1:
					ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
					plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)	
				else: ax.set_xticklabels([])
			
				#---save histograms for overall plot
				chists.append(hist*len(cfilt_inds))
			
				#---plot structure
				ax_struct = axlist_struct[pnum]
				im = plotter2d(ax_struct,mset,dat=mean(mset.surf,axis=0)/10.,
					lognorm=False,cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
					ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
					fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.],
					tickskip=int(round(mset.griddims[0]/6,-1)))
				if pnum < len(params_plot)-1: ax_struct.set_xticklabels([])
				#---height color scale
				axins2 = inset_axes(ax_struct,width="5%",height="100%",loc=3,
					bbox_to_anchor=(1.,0.,1.,1.),
					bbox_transform=ax_struct.transAxes,
					borderpad=0)
				cbar = plt.colorbar(im,cax=axins2,orientation="vertical")
				axins2.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',
					fontsize=fsaxlabel,rotation=270)
				axins2.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
				if pnum == 0 and 0: ax_struct.set_title(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$')
				#---plot hulls according to scheme type
				if re.match('lattice.+',plot_scheme_name):
					abs_pts = mean(nborhoods[cnum],axis=0)[0]
					pt = [abs_pts[i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] for i in range(2)]
					ax_struct.add_patch(plt.Circle((pt[0],pt[1]),
						radius=0.5*mset.griddims[0]/35,color=color,alpha=1.))
				elif plot_scheme_name == 'proteins_controls':
					if nbornames[cnum] != 'oligomer':
						protpts = mean(nborhoods[cnum],axis=0)
						plothull(ax_struct,protpts,mset=mset,subdivide=None,
							c=color,alpha=1.,fill=(False if geog[0] == 'control' else True))
				if plot_scheme_name == 'extrema' and label in ['valley (dynamic)','peak (dynamic)']:
					for protpts in nborhoods[cnum]:
						plothull(ax_struct,protpts,mset=mset,subdivide=None,
							c=color,alpha=1.,fill=(False if geog[0] == 'control' else True),
							radius=0.25*mset.griddims[0]/35)
				if plot_scheme_name == 'extrema' and label in ['valley','peak']:
					pt = [mean(nborhoods[cnum],axis=0)[0][i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] 
						for i in range(2)]
					for shift in [[i,j] for j in [-1,0,1] for i in [-1,0,1]]:
						ptpbc = [pt[i]+shift[i]*mset.griddims[i] for i in range(2)]
						ax_struct.add_patch(plt.Circle((ptpbc[0],ptpbc[1]),radius=2*mset.griddims[0]/35,
							color='w',alpha=1.,lw=2,ec='k'))
						ax_struct.plot([ptpbc[0]],[ptpbc[1]],marker=('_' if label == 'valley' else '+'),
							mec='k',ms=6,mew=2)
				ax_struct.set_xlim((0,mset.griddims[0]-1))
				ax_struct.set_ylim((0,mset.griddims[1]-1))
					
				#---plot extents of curvature
				extdat = sqrt(abs(array([test.data[i][0][4:6] for i in cfilt_inds \
					if type(test.data[i][1]) != list]).flatten()))
				hist,edges = numpy.histogram(extdat,
					range=(0,extent_range),bins=extent_range+1)
				axlist_extent[pnum].plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2)
				if max(hist) > maxpeak_extent: maxpeak_extent = max(hist)
				if pnum == len(params_plot)-1:
					axlist_extent[pnum].get_xaxis().set_major_locator(\
						mpl.ticker.MaxNLocator(prune='both',nbins=4))
					axlist_extent[pnum].set_xlabel('extent '+r'$\mathrm{\sigma_{a,b}\,(nm)}$',
						fontsize=fsaxlabel)
						
				#---computations
				valid_frames = len(hmaxdat)
				result_nums = [
					(round(mean(hmaxdat),4) if valid_frames > 0 else np.nan),
					round(mean(resids),1),
					round(mean(extdat[extdat<2*extent_range]),2),
					round(mean(extdat),2)]
				geog = listlook(test.notes,'geography')
				if geog[0] == 'control':
					casual_name = analysis_descriptors[geog[1]]['label_text']
				else: casual_name = nbornames[cnum]
				fp.write(
					str(casual_name).ljust(12)+\
					str((' ' if result_nums[0] > 0 else '')+str(result_nums[0])).ljust(18)+\
					str('%1.0f'%valid_frames).ljust(8)+\
					str(str(' ' if result_nums[1] > 0 else '')+str(result_nums[1])).ljust(8)+\
					str(str(' ' if result_nums[2] > 0 else '')+str(result_nums[2])).ljust(16)+\
					str(str(' ' if result_nums[3] > 0 else '')+str(result_nums[3])).ljust(25)+'\n')
					
				#---save the mean Hmax if performing the lattice test for later histogram
				if re.match('lattice.+',plot_scheme_name) and not isnan(result_nums[0]):
					lattice_mean_collect[pnum].append(result_nums[0])
	
			#---truncate the legend if too large
			if len(plot_order) < 8: ax.legend(loc='upper left',prop={'size':fsaxlegend_small})
			else:
				which_label = slice(-1,None) if plot_order[0] != 0 else slice(None,1)
				h,l = ax.get_legend_handles_labels()
				ax.legend(h[which_label],l[which_label],loc='upper left',prop={'size':fsaxlegend_small})
			
		
		#---collect the sum of each curve for further analysis
		phists.append(mean(chists,axis=0))
		sysnames = [listlook(params_plot[p][0],'label_text') for p in range(len(params_plot))]
		#---plot settings
		for a in range(len(axlist_extent)):
			ax = axlist[a]
			ax.grid(True)
			ax.set_xlim((-hifilt,hifilt))
			ax.set_ylim((0,1.1*maxpeak))
			ax.set_yticklabels([])
			color,proper_name = dimple_colordict_systems(sysnames[a])
			ax.set_ylabel(proper_name,fontsize=fsaxlabel)
			ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')
		for a in range(len(axlist_extent)):
			ax = axlist_extent[a]
			ax.set_ylim((0,1.1*maxpeak_extent))
			ax.grid(True)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
			ax.set_yticklabels([])
			if a < len(axlist_extent)-1: ax.set_xticklabels([])
			ax.set_xlim((0,extent_range))

	#---clean and save
	print '\nstatus: saving plot'
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	plt.savefig(pickles+'fig-dimple3-'+bigname+\
		'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
		plot_scheme_name+\
		'.png',dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	fp.close()
	print 'status: plot saved'

#---additional summary plot
if 'plot' in routine and re.match('lattice.+',plot_scheme_name):

	#---plot the sum of the distributions for clarity (or weight by quantity of valid frames)
	fig = plt.figure(figsize=(10,10))
	ax1 = plt.subplot(211)
	ax2 = plt.subplot(212)
	nnames = [listlook(params_plot[p][0],'label_text') for p in range(len(phists))]
	for p in range(len(phists)):
		color,proper_name = dimple_colordict_systems(nnames[p])
		status('curve = '+str(p+1)+'/'+str(len(params_plot)))
		ax1.plot(hmax_mids,phists[p],lw=2,color=color,
			label=proper_name+' ('+str(round(sum(phists[p]*hmax_mids)/sum(phists[p]),3))+\
			r'$\mathsf{\,(nm^{-1})}$'+')')
		hist,edges = numpy.histogram(lattice_mean_collect[p],
			range=(-hifilt-hist_step/2,hifilt+hist_step/2),bins=hmax_nbins)
		ax2.bar(1./2*(edges[1:]+edges[:-1]),hist*sum(phists[p])/sum(hist),'-',c=color,lw=2,
			label=proper_name+' ('+str(round(mean(lattice_mean_collect[p]),3))+r'$\mathsf{\,nm^{-1}}$'+')')
	#---plot parameters
	for ax in [ax1,ax2]:
		ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')	
		ax.legend(prop={'size':fsaxlegend_small})
		ax.grid(True)
		ax.set_xlim((-hifilt,hifilt))
		ax.set_yticklabels([])
		ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')
		ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)	
	ax1.set_ylabel('summed distributions',fontsize=fsaxlabel)
	ax2.set_ylabel('distribution of means',fontsize=fsaxlabel)
	ax1.set_title('lattice test',fontsize=fsaxtitle)
	
	#---clean and save
	print '\nstatus: saving plot'
	plt.savefig(pickles+'fig-dimple3-'+bigname+\
		'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
		plot_scheme_name+'.summary'\
		'.png',dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	print 'status: plot saved'
	

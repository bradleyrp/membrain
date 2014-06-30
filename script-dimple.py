#!/usr/bin/python -i

interact = True
from membrainrunner import *
execfile('locations.py')

#---imports
from scipy.optimize import leastsq
import matplotlib.image as mpimg

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

#---analysis plan
analysis_names = [
	'v616-210000-310000-200',
	'v614-120000-220000-200',
	'v612-75000-175000-200',
	'v550-400000-500000-160',
	'v614-40000-140000-200',
	'v612-10000-80000-200',
	'v550-300000-400000-200',
	'v700-500000-600000-200',
	'v701-60000-160000-200',
	][1:4]
routine = [
	'compute_dimple',
	'plot_dimple',
	'plot_dimple_pub',
	'plot_radar',
	'compute_topog',
	'plot_topog',
	][2:3]
bigname = 'v616-v614-v612-v550-ver2'
bigname = 'v614-v612-v550-IET-revise'

#---available dimple plots
dimple_plot_type = [
	'lattice6',
	'extrema',
	'proteins_controls',
	'proteins_controls_pub',
	'extrema_meta',
	][3]
	
#---available radar plots
radar_plot_type = ['radar_proteins','radar_extrema'][-1]

#---available topography plots
topog_plot_type = [
	'extrema',
	][0]
	
show_plots = False
	
#---DOWNSTREAM
#-------------------------------------------------------------------------------------------------------------
	
#---set downstream parameters
if sum([i in routine for i in ['compute_dimple','compute_topog']]) > 1:
	raise Exception('except: contradicting computations')
#---settings for test types
if any([i in routine for i in ['compute_dimple','plot_dimple','plot_dimple_pub']]): 
	execfile('script-dimple-specs.py')
	pklprefix = 'pkl.dimple3.'
	testtype = 'dimple'
elif any([i in routine for i in ['compute_topog','plot_topog']]): 
	execfile('specs-topography.py')
	pklprefix = 'pkl.topog3.'
	testtype = 'topography'
	
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams.update({'font.size': 14})

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---NOTE: the function gauss2dh returns a value equal to 2H, which is corrected with the following factor
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

def dimple_save_data(aname):
	'''Saves the global dimple3 result data object for a particular system name (name msdats[anum]).'''
	for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
	anum = analysis_names.index(aname)
	filespec = specname_guess(sysname,trajsel)
	pklname = pklprefix+spacetag+filespec+'.pkl'
	pickledump(msdats[anum],pklname,directory=pickles)

def listlook(data,name):
	'''Treat a list like a dictionary and look up a value.'''
	for i in data:
		if i[0] == name: return i[1]

def dimple_test_lookup(notes,aname):
	'''Find the indices corresponding to a particular test.'''
	anum = analysis_names.index(aname)
	inds = [j for j in range(len(msdats[anum])) if all([i in msdats[anum][j].notes for i in notes]) == True]
	#inds = [j for j in range(len(msdats[anum])) if all([i in notes for i in msdats[anum][j].notes]) == True]
	if inds != []: return inds
	else: return False

def prepare_experiments():	
	'''Prepare a list of experiments according to any requested sweeps.'''
	list_params_expts = []
	keylist = params_sweeps.keys()
	sweeps = [params_sweeps[key] for key in keylist]
	param_combos = list(itertools.product(*sweeps))
	for new_expt in param_combos:
		expt_copy = params_expt.copy()
		for key in keylist:
			expt_copy[key] = new_expt[keylist.index(key)]
		list_params_expts.append(expt_copy)
	return list_params_expts
	
def compile_expt_notes(aname,expt,analysis_descriptors):
	'''Compile a "notes" object with all parameters for a test+system combination.'''
	notes_list = []
	for i in analysis_descriptors[aname]:
		notes_list.append([i,(analysis_descriptors[aname])[i]])
	for key in expt.keys():
		notes_list.append([key,expt[key]])
	notes_list.append(['callsign',aname])
	return notes_list
	
def compile_plot_params(plotname,analysis_descriptors):
	'''Selects a plotname from the master list and compiles experiment/system combinations for plotting.'''
	params_plot = []
	for i in params_plot_master[params_plot_master_names.index(plotname)]:
		subplot = []
		for j in i:
			subplot.append(compile_expt_notes(j[0],j[1],analysis_descriptors))
		params_plot.append(subplot)
	return params_plot

def dimple_generate_neighborhood(notes):
	'''Generate neighborhood points.'''
	request_geography = listlook(notes,'geography')
	#---select the surface
	aname = listlook(notes,'callsign')
	anum = analysis_names.index(aname)
	mset = msets[anum]

	#---select frames
	framecounts = listlook(notes,'framecounts')
	if framecounts != None and type(framecounts) != slice: frame_select = range(framecounts)
	elif type(framecounts) == slice: frame_select = range(len(mset.surf))[framecounts]
	else: frame_select = range(len(mset.surf))


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
		
def dimple_colordict(ask,listing=None,scheme=None,nprots=None,reserve=True):
	'''Color lookup for plotting protein hulls.'''
	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in [0,1,2,3,4,6,7,8]]
	clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in [0,1,2,3,4,6,8,7]]
	clrsdict_specific = {'red':clrs[0],'blue':clrs[1],'black':'k'}
	clrsdict_leftovers = list(array(clrs)[range(2,len(clrs))])
	if scheme != None and re.match('lattice.+',scheme):
		clrsdict_leftovers = [brewer2mpl.get_map('Set3','qualitative',12).mpl_colors[i] for i in range(12)]
	if reserve: reserved_colors = ['max','min,','max (dynamic)','min (dynamic)','oligomer']
	else: 
		reserved_colors = []
		clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] for i in [0,1,2,3,4,6,7,8]]
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
	
def positive_negative_areas(mset,neighbor_defn=None):
	'''Function which gets all surface points within a buffer of another set of points.'''
	posneg = []
	if neighbor_defn == None: slicelen = len(mset.protein[fr])
	else: slicelen = len(neighbor_defn)
	slicelen = min([slicelen,len(mset.surf)])
	for fr in range(slicelen):
		status('status: computing areas ',i=fr,looplen=len(mset.surf))
		ptsa = unzipgrid(mset.surf[fr],vecs=mset.vec(fr))
		if neighbor_defn == None: ptsb = mset.protein[fr]
		else: ptsb = neighbor_defn[fr]
		dmat = scipy.spatial.distance.cdist(ptsa[:,:2],ptsb[:,:2])
		nhood = ptsa[where(dmat.min(axis=1)<100)[0]]
		posneg.append([sum(nhood[:,2]>0),sum(nhood[:,2]<0)])
	return array(posneg)
	
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
		
#---NEEDS COMMENT
if 'posnegs' not in globals():
	posneg_select = ['','-4area'][1]
	#--- ???????????????????????/ NEEDS COMMENTS !!!
	posnegs = []
	#---explicit requirement to use the 4xENTH neighborhood for all systems
	for m in range(3):
		mset = msets[m]
		if posneg_select == '':
			if m == 0: neighbor_defn = msets[analysis_names.index('v614-120000-220000-200')].protein
			elif m > 0: neighbor_defn = msets[analysis_names.index('v612-75000-175000-200')].protein
		elif posneg_select == '-4area':
			neighbor_defn = msets[analysis_names.index('v614-120000-220000-200')].protein
		posneg = positive_negative_areas(mset,neighbor_defn=neighbor_defn)
		posnegs.append(posneg)
		
#---compute fits
if 'compute_dimple' in routine or 'compute_topog' in routine:
	
	#---load the current set of msdats which may already have some of the requested tests in the sweep
	msdats = [[] for i in range(len(analysis_names))]
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		anum = analysis_names.index(aname)
		#---load the pickle to see if it exists	
		filespec = specname_guess(sysname,trajsel)
		pklname = pklprefix+spacetag+filespec+'.pkl'
		msdat = unpickle(pickles+pklname)
		if msdat == None: msdats[anum] = []
		else: msdats[anum] = msdat
	#---prepare a list of experiments according to sweeps
	list_params_expts = prepare_experiments()

	#---loop over experiments and systems and combine valid experiments into a list of notes objects
	list_expts = []
	for expt in list_params_expts:
		for aname in analysis_names:
			notes = compile_expt_notes(aname,expt,analysis_descriptors)
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
		elif type(framecounts) == slice: frame_select = range(len(mset.surf))[framecounts]
		else: frame_select = range(len(mset.surf))

		#---loop over neighborhoods within each test type
		for nborhood_ind in range(len(nborhoods)):
			nborhood = nborhoods[nborhood_ind]
		
			#---check to see if this calculation was already completed
			notes_nbor = notes+[['neighborhood_name',nbornames[nborhood_ind]]]
			needs_test = (False == dimple_test_lookup(notes_nbor,aname))

			#---unpack the notes
			for i in notes_nbor: vars()[i[0]] = i[1]

			#---dimple fitting calculation
			if needs_test and testtype == 'dimple':
			
				#---prepare the result file
				result_data = MembraneData('dimple3',label=sysname)
				result_data.notes = notes_nbor
					
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

			#---topography calculation
			if needs_test and testtype == 'topography':

				#---prepare the result file
				result_data = MembraneData('dimple3',label=sysname)
				result_data.notes = notes_nbor	

				#---catalog bilayer heights in the neighborhood
				heights = []
				for fr in frame_select:
					status(str(fr))
					status('neighborhood = '+str(nborhood_ind+1)+'/'+str(len(nborhoods))+' frame = '+str(fr))
					#---replicate surface points under periodic boundary conditions
					surfpts = mset.wrappbc(mset.unzipgrid(mset.surf[fr],vecs=mset.vec(0)),
						vecs=mset.vec(fr),mode='nine')
					#---find minimum distances
					cd = scipy.spatial.distance.cdist(nborhood[fr],surfpts[:,:2])
					tmp = array([np.min(cd,axis=0),surfpts[:,2]]).T
					selected = where(tmp[:,0]<cutoff)[0]
					target = surfpts[selected]
					heights.append(list(target[:,2]))

				#---populate the result object with the dimple fit results
				for i in range(len(heights)):
					result_data.add([heights[i]],[frame_select[i]])
				msdats[anum].append(result_data)
				dimple_save_data(aname)
				

#---MAIN, PLOT
#-------------------------------------------------------------------------------------------------------------

#---plot summary
if 'plot_dimple' in routine:

	if 'msdats' not in globals():
		#---load the current set of msdats
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---load the pickle to see if it exists	
			filespec = specname_guess(sysname,trajsel)
			pklname = pklprefix+spacetag+filespec+'.pkl'
			msdat = unpickle(pickles+pklname)
			if msdat == None: msdats[anum] = []
			else: msdats[anum] = msdat

	'''
	This plotting procedure is the de-facto analysis of the dimple results, and also produces summary tables.
	I attempted to move it to an in-place function.
	Since I unpack analysis_descriptors and params_plot into global variables, I tried updating locals().
	Note that this is not preferred. You can pass a dictionary of variables to a local function like this:
		def funct(a,b): print a
		somedict = {'a':2,'b':3}
		funct(**somedict)
	This doesn't work because unlike globals(), locals() is not a true dictionary, due to optimizations.
	So I skipped putting the plot function in-place and added an external loop for meta-analysis.
	'''
	
	#---meta-analysis loop
	params_plot = compile_plot_params(dimple_plot_type,analysis_descriptors)
	if type(params_plot_settings) == dict: params_plot_settings_loop = [params_plot_settings]
	else: params_plot_settings_loop = params_plot_settings
	#---if the meta-analysis was already complete, skip ahead
	if dimple_plot_type == 'extrema_meta':
		dimple_meta = unpickle(pickles+'fig-dimple3-meta-extrema-'+bigname+'.pkl')
		if dimple_meta != None: params_plot_settings_loop = []
		else: dimple_meta = []
	for params_plot_settings in params_plot_settings_loop:
	
		#---only do meta analysis if the params_plot is multidimensional
		if len(shape(params_plot)) > 1:
			dimple_meta_set = [[[] for i in range(shape(params_plot)[1])] 
			for j in range(shape(params_plot)[0])]
		
		#---compile the plot specifications
		for i in params_plot_settings: globals()[i] = params_plot_settings[i]
	
		#---prepare plot panels	
		fig = plt.figure(figsize=(10,len(params_plot)*2))
		gs = gridspec.GridSpec(len(params_plot),3,wspace=0.0,hspace=0.0)
		gs.update(left=0.0,right=0.7)
		gs2 = gridspec.GridSpec(len(params_plot),1,wspace=0.0,hspace=0.0)
		gs2.update(left=0.75,right=1.0)
		axlist,axlist_struct,axlist_extent = [],[],[]
		fp = open(pickles+'calc-dimple3-'+bigname+\
			'.cut'+str(listlook(params_plot[0][0],'cutoff'))+spacetag+\
			'filter-'+filter_type+'.'+\
			dimple_plot_type+'.txt','w')

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
				for i in analysis_descriptors[aname]: globals()[i] = (analysis_descriptors[aname])[i]
				nborhoods,nbornames = dimple_generate_neighborhood(expt)
			
				#---write header for the calculations text output
				if dimple_plot_type != 'extrema' or gnum == 0:
					fp.write('\nsystem = '+str(aname)+'\t'+(analysis_descriptors[aname])['label_text']+'\n\n')
					table_top = 'name                '+'<H_max> (nm^-1)   '+'<H_max> scaled    '+\
						'frames  '+'RMSD(A) '+'sigma_a,b (nm)  '+'sigma_a,b unfiltered (nm)'
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
						if listlook(test.notes,'nprots') == 1 and dimple_plot_type != 'extrema' and \
							re.match('lattice.+',dimple_plot_type) == None:
							nprots = 1
						else: nprots = None
						nprots = (listlook(test.notes,'nprots') if \
								re.match('lattice.+',dimple_plot_type) == None else None)
						color = dimple_colordict(nnames[cnum],listing=nnames,
							scheme=dimple_plot_type,nprots=nprots)
					if label in ['peak (dynamic)','valley (dynamic)']: alpha = 0.5
					else: alpha = 1.
					status('panel = '+str(pnum+1)+'/'+str(len(params_plot))+
						' curve = '+str(cnum+1)+'/'+str(len(dat_inds)))
					#---standard filtering
					if filter_type == 'std':
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
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					#---alternate filtering 1
					#---note that the "mod1" handle will never change
					#---this filters for magnitudes within smallfilt and hifilt and dimples within cutoff 
					#---...of the COM of the full neighborhood so that for point neighborhoods with buffer
					#---...this method only counts dimples that are inside of the neighborhood
					elif filter_type == 'mod1':
						vecs = mean(mset.vecs,axis=0)
						fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
						#---dropped the following test for being the box because PBCs were used on the 
						if 0: inbox_inds = array([all([i[0][2] < vecs[0] and i[0][3] < vecs[1] and 
							i[0][2] > 0. and i[0][3] > 0.]) for i in test.data])
						#---hack only works for pointsize neighborhoods
						inbox_inds = array([linalg.norm(array([i[0][2],i[0][2]])-\
							mean(i[-1],axis=0)[:2])/10.<listlook(test.notes,'cutoff') for i in test.data])
						hmaxraw = []
						for i in range(len(test.data)):
							if fitted_inds[i]:
								z0,c0,x0,y0,sx,sy,th = test.data[i][0]
								hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
							else: hmaxraw.append(0)
						magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
							for i in range(len(hmaxraw))])
						cfilt_inds = where(array(1*magfilter+1*inbox_inds+1*fitted_inds)==3)[0]
						hmaxdat = array(hmaxraw)[cfilt_inds]
						params = test.get(['type','params'])[cfilt_inds]
						target_zones = test.get(['type','target_zones'])[cfilt_inds]
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					#---alternate filtering 2
					#---note that the "mod2" handle will never change
					#---filters for magnitudes within smallfilt and hifilt and dimples in the neighborhood
					elif filter_type == 'mod2':
						vecs = mean(mset.vecs,axis=0)
						fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
						target_zones = test.get(['type','target_zones'])
						hmaxraw = []
						center_near_nborhood = [False for i in range(len(test.data))]
						for i in range(len(test.data)):
							if fitted_inds[i]:
								z0,c0,x0,y0,sx,sy,th = test.data[i][0]
								hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
								center_near_nborhood[i] = scipy.spatial.distance.cdist(target_zones[i][:,:2],
									[[x0,y0]]).min() < sqrt(2)*10.
							else: hmaxraw.append(0)
						center_near_nborhood = array(center_near_nborhood)
						magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
							for i in range(len(hmaxraw))])
						cfilt_inds = where(array(1*magfilter+1*center_near_nborhood+1*fitted_inds)==3)[0]
						hmaxdat = array(hmaxraw)[cfilt_inds]
						params = test.get(['type','params'])[cfilt_inds]
						target_zones = test.get(['type','target_zones'])[cfilt_inds]
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					else: print 'error: filter modification is underspecified'
					'''
					NOTE THAT I CHECKED THE MOD 2 METHOD WITH THE FOLLOWING CODEBLOCK
					if 0:
						params = test.get(['type','params'])
						target_zones = test.get(['type','target_zones'])
						thelist = list(where(center_near_nborhood==False)[0])
						for i in thelist:
							x0,y0 = params[i][2:4]
							dat = scipy.spatial.distance.cdist(target_zones[i][:,:2],[[x0,y0]])
							j = argmin(dat[:,0])
							#j = list(dat[:,0]).index(sort(dat[:,0])[0])
							meshpoints(array(list(params[i][2:4])+[0.0]),scale_factor=20)
							meshpoints(target_zones[i],scale_factor=10,color=(1,1,1))
							meshpoints(target_zones[i][j],scale_factor=10,color=(1,0,1))
					'''
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
					if re.match('lattice.+',dimple_plot_type):
						abs_pts = mean(nborhoods[cnum],axis=0)[0]
						pt = [abs_pts[i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] for i in range(2)]
						ax_struct.add_patch(plt.Circle((pt[0],pt[1]),
							radius=0.5*mset.griddims[0]/35,color=color,alpha=1.))
					elif dimple_plot_type == 'proteins_controls':
						if nbornames[cnum] != 'oligomer':
							protpts = mean(nborhoods[cnum],axis=0)
							plothull(ax_struct,protpts,mset=mset,subdivide=None,
								c=color,alpha=1.,fill=(False if geog[0] == 'control' else True))
					if dimple_plot_type == 'extrema' and label in ['valley (dynamic)','peak (dynamic)']:
						for protpts in nborhoods[cnum]:
							plothull(ax_struct,protpts,mset=mset,subdivide=None,
								c=color,alpha=1.,fill=(False if geog[0] == 'control' else True),
								radius=0.25*mset.griddims[0]/35)
					if dimple_plot_type == 'extrema' and label in ['valley','peak']:
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
					scaled_mean = float(valid_frames)/len(test.data)*mean(hmaxdat)
					result_nums = [
						(round(mean(hmaxdat),4) if valid_frames > 0 else np.nan),
						(round(scaled_mean,4) if valid_frames > 0 else np.nan),
						round(mean(resids),1),
						round(mean(extdat[extdat<2*extent_range]),2),
						round(mean(extdat),2),
						valid_frames]
					print result_nums
					geog = listlook(test.notes,'geography')
					if geog[0] == 'control':
						casual_name = analysis_descriptors[geog[1]]['label_text']
					else: casual_name = nbornames[cnum]
					fp.write(
						str(casual_name).ljust(20)+\
						str((' ' if result_nums[0] > 0 else '')+str(result_nums[0])).ljust(18)+\
						str((' ' if result_nums[1] > 0 else '')+str(result_nums[1])).ljust(18)+\
						str('%1.0f'%valid_frames).ljust(8)+\
						str(str(' ' if result_nums[2] > 0 else '')+str(result_nums[2])).ljust(8)+\
						str(str(' ' if result_nums[3] > 0 else '')+str(result_nums[3])).ljust(16)+\
						str(str(' ' if result_nums[4] > 0 else '')+str(result_nums[4])).ljust(25)+'\n')
					#---save the mean Hmax if performing the lattice test for later histogram
					if re.match('lattice.+',dimple_plot_type) and not isnan(result_nums[0]):
						lattice_mean_collect[pnum].append(result_nums[0])
						
					#--save for meta-analysis
					#---only do meta analysis if the params_plot is multidimensional
					if len(shape(params_plot)) > 1:
						dimple_meta_set[anum][gnum] = result_nums

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
		if dimple_plot_type != 'extrema_meta':
			print '\nstatus: saving plot'
			fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
			plt.savefig(pickles+'fig-dimple3-'+bigname+\
				'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
				'filter-'+filter_type+'.'+\
				dimple_plot_type+\
				'.png',dpi=300,bbox_inches='tight')
			if show_plots: plt.show()
			plt.close()
			fp.close()
			print 'status: plot saved'
		
		#---only do meta analysis if the params_plot is multidimensional
		if len(shape(params_plot)) > 1:
			dimple_meta.append(dimple_meta_set)
		
	#---dump meta data
	if dimple_plot_type == 'extrema_meta':
		pickledump(dimple_meta,'fig-dimple3-meta-extrema-'+bigname+'.pkl',directory=pickles)
		
#---PUBLICATION PLOT
#-------------------------------------------------------------------------------------------------------------

'''Revised from the main plot above for IET paper.'''

#---plot summary
if 'plot_dimple_pub' in routine:

	'''
	This is similar to the plot_dimple algorithm with some modifications for the IET paper.
	pubstyle2 puts the oligomer behind the monomers if there are any
	it also removes the fill from the monomers so they don't block it
	and it also changes the color scheme to a random one
	'''
	pubstyle2 = True
	showstruct = True
	#---set to false to use all available colors (set to True to avoid blue/red for valley/peak
	color_override = True

	if 'msdats' not in globals():
		#---load the current set of msdats
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---load the pickle to see if it exists	
			filespec = specname_guess(sysname,trajsel)
			pklname = pklprefix+spacetag+filespec+'.pkl'
			msdat = unpickle(pickles+pklname)
			if msdat == None: msdats[anum] = []
			else: msdats[anum] = msdat

	#---meta-analysis loop
	params_plot = compile_plot_params(dimple_plot_type,analysis_descriptors)
	if type(params_plot_settings) == dict: params_plot_settings_loop = [params_plot_settings]
	else: params_plot_settings_loop = params_plot_settings
	#---if the meta-analysis was already complete, skip ahead
	if dimple_plot_type == 'extrema_meta':
		dimple_meta = unpickle(pickles+'fig-dimple3-meta-extrema-'+bigname+'.pkl')
		if dimple_meta != None: params_plot_settings_loop = []
		else: dimple_meta = []
	for params_plot_settings in params_plot_settings_loop:
		
		#---only do meta analysis if the params_plot is multidimensional
		if len(shape(params_plot)) > 1:
			dimple_meta_set = [[[] for i in range(shape(params_plot)[1])] 
			for j in range(shape(params_plot)[0])]
		
		#---compile the plot specifications
		for i in params_plot_settings: globals()[i] = params_plot_settings[i]
	
		#---prepare plot panels	
		fig = plt.figure(figsize=(8,len(params_plot)*2))
		gs = gridspec.GridSpec(len(params_plot),3,wspace=0.0,hspace=0.0)
		gs.update(left=0.0,right=0.67)
		gs2 = gridspec.GridSpec(len(params_plot),3,wspace=0.0,hspace=0.0)
		gs2.update(left=0.69,right=1.0)
		#axlist,axlist_struct,axlist_extent = [],[],[]
		axlist,axlist_areas,axlist_areas2,axlist_extent = [],[],[],[]
		fp = open(pickles+'calc-dimple3-'+bigname+\
			'.cut'+str(listlook(params_plot[0][0],'cutoff'))+spacetag+\
			'filter-'+filter_type+'.'+\
			dimple_plot_type+'.txt','w')

		#---global plot specifications
		extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())])
			for mset in msets])
		maxpeak,maxpeak_extent = 0,0
		maxpeak_areas,maxfreq = 0,0
		phists = []
		lattice_mean_collect = [[] for j in range(len(params_plot))]
		
		#---loop over rows
		for pnum in range(len(params_plot)):
			ax = plt.subplot(gs[pnum,0:2])
			axlist.append(ax)
			#ax_struct = plt.subplot(gs2[pnum])
			#axlist_struct.append(ax_struct)
			ax_areas = plt.subplot(gs2[pnum,0:2])
			ax_areas2 = plt.subplot(gs2[pnum,2])
			axlist_areas.append(ax_areas)
			axlist_areas2.append(ax_areas2)
			ax_extent = plt.subplot(gs[pnum,2])
			axlist_extent.append(ax_extent)
			nnames_control = []
			
			#---area plots, outside of the loop over plots-within-a-panel, because only one can be plotted
			ax = axlist_areas[pnum]
			#---codeblock from script-reproc-dimple-backup
			poz_area_counts = posnegs[pnum][:,0]
			neg_area_counts = posnegs[pnum][:,1]
			area_per_tile = product(mean(msets[pnum].vecs,axis=0)[:2])/100./\
				((msets[pnum].griddims[0]-1)*(msets[pnum].griddims[1]-1))
			posarea = array([area_per_tile*poz_area_counts[i] for i in range(len(poz_area_counts))])
			negarea = array([area_per_tile*neg_area_counts[i] for i in range(len(neg_area_counts))])
			ax.plot(posarea,'r-',label='$z>0$',lw=1,alpha=0.8)
			ax.plot(negarea,'b-',label='$z<0$',lw=1,alpha=0.8)
			t = range(len(poz_area_counts))
			ax.fill_between(t, negarea,posarea, facecolor='b',alpha=0.25,where=negarea>posarea)
			ax.fill_between(t, posarea,negarea, facecolor='r',alpha=0.25,where=negarea<posarea)
			if max(max(posarea),max(negarea)) > maxpeak_areas: maxpeak_areas = max(max(posarea),max(negarea))

			nbins_areas = 40
			minval_areas = 0
			maxval_areas = maxpeak_areas
			ax = axlist_areas2[pnum]
			hist0,binedge0 = numpy.histogram(posarea,bins=nbins_areas,normed=False,
				weights=[1./len(negarea) for i in negarea],
				range=(minval_areas,maxval_areas*1.35))
			mid0 = (binedge0[1:]+binedge0[:-1])/2
			ax.plot(hist0,mid0,c='r',alpha=1.,lw=1,label='$z>0$')
			hist1,binedge1 = numpy.histogram(negarea,bins=nbins_areas,normed=False,
				weights=[1./len(negarea) for i in negarea],
				range=(minval_areas,maxval_areas*1.35))
			mid1 = (binedge1[1:]+binedge1[:-1])/2
			ax.plot(hist1,mid1,c='b',alpha=1.,lw=1,label='$z<0$')
			ax.fill_betweenx(mid0,[0 for i in mid0],hist0,facecolor='r',
				alpha=0.2)
			ax.fill_betweenx(mid1,[0 for i in mid1],hist1,facecolor='b',
				alpha=0.2)
			if max(max(hist0),max(hist1)) > maxfreq: maxfreq = max(max(hist0),max(hist1))
			
			#---loop over plots within a panel
			chists = []
			for gnum in range(len(params_plot[pnum])):
				expt = params_plot[pnum][gnum]
				aname = listlook(expt,'callsign')
				anum = analysis_names.index(aname)
				mset = msets[anum]
				for i in analysis_descriptors[aname]: globals()[i] = (analysis_descriptors[aname])[i]
				nborhoods,nbornames = dimple_generate_neighborhood(expt)
			
				#---write header for the calculations text output
				if dimple_plot_type != 'extrema' or gnum == 0:
					fp.write('\nsystem = '+str(aname)+'\t'+(analysis_descriptors[aname])
						['label_text']+'\n\n')
					table_top = 'name'.ljust(20)+\
						'<H_max> (nm^-1)'.ljust(18)+\
						'<H_max> scaled'.ljust(18)+\
						'frames'.ljust(8)+'RMSD(A)'.ljust(10)+\
						'sigma_a'.ljust(10)+'sigma_b'.ljust(10)+\
						'sigma_a,b (nm)'.ljust(16)+\
						'sigma_a,b unfiltered (nm)'.ljust(25)
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
				if 'oligomer' in nbornames and pubstyle2: plot_order = [0]+range(len(dat_inds))[1:]
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
						color = 'k'
					else: 
						label = listlook(test.notes,'neighborhood_name')
						if listlook(test.notes,'nprots') == 1 and dimple_plot_type != 'extrema' and \
							re.match('lattice.+',dimple_plot_type) == None:
							nprots = 1
						else: nprots = None
						nprots = (listlook(test.notes,'nprots') if \
								re.match('lattice.+',dimple_plot_type) == None else None)
						color = dimple_colordict(nnames[cnum],listing=nnames,
							scheme=dimple_plot_type,nprots=nprots)
					#--- !!!
					if color_override and len(plot_order)>1 and nbornames[cnum] != 'oligomer':
						color = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[i] 
							for i in [0,1,2,3,4,6,7,8]][cnum]
					
					if label in ['peak (dynamic)','valley (dynamic)']: alpha = 0.5
					else: alpha = 1.
					status('panel = '+str(pnum+1)+'/'+str(len(params_plot))+
						' curve = '+str(cnum+1)+'/'+str(len(dat_inds)))
					#---standard filtering
					if filter_type == 'std':
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
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					#---alternate filtering 1
					#---note that the "mod1" handle will never change
					#---this filters for magnitudes within smallfilt and hifilt and dimples within cutoff 
					#---...of the COM of the full neighborhood so that for point neighborhoods with buffer
					#---...this method only counts dimples that are inside of the neighborhood
					elif filter_type == 'mod1':
						vecs = mean(mset.vecs,axis=0)
						fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
						#---dropped the following test for being the box because PBCs were used on the 
						if 0: inbox_inds = array([all([i[0][2] < vecs[0] and i[0][3] < vecs[1] and 
							i[0][2] > 0. and i[0][3] > 0.]) for i in test.data])
						#---hack only works for pointsize neighborhoods
						inbox_inds = array([linalg.norm(array([i[0][2],i[0][2]])-\
							mean(i[-1],axis=0)[:2])/10.<listlook(test.notes,'cutoff') for i in test.data])
						hmaxraw = []
						for i in range(len(test.data)):
							if fitted_inds[i]:
								z0,c0,x0,y0,sx,sy,th = test.data[i][0]
								hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
							else: hmaxraw.append(0)
						magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
							for i in range(len(hmaxraw))])
						cfilt_inds = where(array(1*magfilter+1*inbox_inds+1*fitted_inds)==3)[0]
						hmaxdat = array(hmaxraw)[cfilt_inds]
						params = test.get(['type','params'])[cfilt_inds]
						target_zones = test.get(['type','target_zones'])[cfilt_inds]
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					#---alternate filtering 2
					#---note that the "mod2" handle will never change
					#---filters for magnitudes within smallfilt and hifilt and dimples in the neighborhood
					elif filter_type == 'mod2':
						vecs = mean(mset.vecs,axis=0)
						fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
						target_zones = test.get(['type','target_zones'])
						hmaxraw = []
						center_near_nborhood = [False for i in range(len(test.data))]
						for i in range(len(test.data)):
							if fitted_inds[i]:
								z0,c0,x0,y0,sx,sy,th = test.data[i][0]
								hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
								center_near_nborhood[i] = scipy.spatial.distance.cdist(target_zones[i][:,:2],
									[[x0,y0]]).min() < sqrt(2)*10.
							else: hmaxraw.append(0)
						center_near_nborhood = array(center_near_nborhood)
						magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
							for i in range(len(hmaxraw))])
						cfilt_inds = where(array(1*magfilter+1*center_near_nborhood+1*fitted_inds)==3)[0]
						hmaxdat = array(hmaxraw)[cfilt_inds]
						params = test.get(['type','params'])[cfilt_inds]
						target_zones = test.get(['type','target_zones'])[cfilt_inds]
						resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
							for i in target_zones[j]])) for j in range(len(target_zones))]
					else: print 'error: filter modification is underspecified'
					'''
					NOTE THAT I CHECKED THE MOD 2 METHOD WITH THE FOLLOWING CODEBLOCK
					if 0:
						params = test.get(['type','params'])
						target_zones = test.get(['type','target_zones'])
						thelist = list(where(center_near_nborhood==False)[0])
						for i in thelist:
							x0,y0 = params[i][2:4]
							dat = scipy.spatial.distance.cdist(target_zones[i][:,:2],[[x0,y0]])
							j = argmin(dat[:,0])
							#j = list(dat[:,0]).index(sort(dat[:,0])[0])
							meshpoints(array(list(params[i][2:4])+[0.0]),scale_factor=20)
							meshpoints(target_zones[i],scale_factor=10,color=(1,1,1))
							meshpoints(target_zones[i][j],scale_factor=10,color=(1,0,1))
					'''
					#---filtered data
					hist,edges = numpy.histogram(hmaxdat[abs(hmaxdat)>smallfilt],
						range=(-hifilt-hist_step/2,hifilt+hist_step/2),bins=hmax_nbins)
					if max(hist) > maxpeak: maxpeak = max(hist)
					#ax.plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2,label=label,alpha=alpha)
					if pubstyle2 and nbornames[cnum] != 'oligomer' and len(nbornames) > 1:
						ax.plot(1./2*(edges[1:]+edges[:-1]),hist,drawstyle='steps-mid',color=color,
							lw=2,label=label,alpha=1.)
					else:
						ax.bar(edges[:-1],hist,width=edges[1]-edges[0],color=color,
							lw=0,label=label,alpha=1.0)
					hmax_mids = 1./2*(edges[1:]+edges[:-1])
					if pnum == len(params_plot)-1:
						ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
						plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)	

						#ax.set_xticks(arange(maxhrange[0],maxhrange[1]+0.001,maxhstep))
						ax.spines['bottom'].set_position(('outward', 10))
						ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
						second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
						ax.spines['second_bottom'] = second_bottom
					else: ax.set_xticklabels([])
					
					#---modified for IET revisions stage 2014.06.13
					letterlist = ['a','b','c','d','e','f','g','h','i','j','k','l'][0:]
					ax.text(0.9,0.86,r'$\mathbf{('+letterlist[pnum].capitalize()+')}$',
						transform=ax.transAxes,fontsize=14)
	
					#---save histograms for overall plot
					chists.append(hist*len(cfilt_inds))
					'''
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
					if pnum == 0 and 0: ax_struct.set_title(
						r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$')
					#---plot hulls according to scheme type
					if re.match('lattice.+',dimple_plot_type):
						abs_pts = mean(nborhoods[cnum],axis=0)[0]
						pt = [abs_pts[i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] for i in range(2)]
						ax_struct.add_patch(plt.Circle((pt[0],pt[1]),
							radius=0.5*mset.griddims[0]/35,color=color,alpha=1.))
					elif dimple_plot_type == 'proteins_controls' or \
						dimple_plot_type == 'proteins_controls_pub':
						if nbornames[cnum] != 'oligomer':
							protpts = mean(nborhoods[cnum],axis=0)
							plothull(ax_struct,protpts,mset=mset,subdivide=None,
								c=color,alpha=1.,fill=(True if geog[0] == 'control' else True))
					if dimple_plot_type == 'extrema' and label in ['valley (dynamic)','peak (dynamic)']:
						for protpts in nborhoods[cnum]:
							plothull(ax_struct,protpts,mset=mset,subdivide=None,
								c=color,alpha=1.,fill=(True if geog[0] == 'control' else True),
								radius=0.25*mset.griddims[0]/35)
					if dimple_plot_type == 'extrema' and label in ['valley','peak']:
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
					'''
					
					#---plot extents of curvature
					extdat = sqrt(abs(array([test.data[i][0][4:6] for i in cfilt_inds 
						if type(test.data[i][1]) != list]).flatten()))
					sigmas0 = array([test.data[i][0][4:6] for i in cfilt_inds 
						if type(test.data[i][1]) != list])[:,0]
					sigmas1 = array([test.data[i][0][4:6] for i in cfilt_inds 
						if type(test.data[i][1]) != list])[:,1]
					sigmas = array([test.data[i][0][4:6] for i in cfilt_inds 
						if type(test.data[i][1]) != list])
					hist,edges = numpy.histogram(extdat,
						range=(0,extent_range),bins=extent_range+1)
					#axlist_extent[pnum].plot(1./2*(edges[1:]+edges[:-1]),hist,'-',c=color,lw=2)
					#axlist_extent[pnum].bar(edges[:-1],hist,width=edges[1]-edges[0],color=color,lw=0)
					ax = axlist_extent[pnum]
					if pubstyle2 and nbornames[cnum] != 'oligomer' and len(nbornames) > 1:
						ax.plot(1./2*(edges[1:]+edges[:-1]),hist,drawstyle='steps-mid',color=color,
							lw=2,label=label,alpha=1.)
					else:
						ax.bar(edges[:-1],hist,width=edges[1]-edges[0],color=color,
							lw=0,label=label,alpha=1.0)				
					if max(hist) > maxpeak_extent: maxpeak_extent = max(hist)
					if pnum == len(params_plot)-1:
						axlist_extent[pnum].get_xaxis().set_major_locator(\
							mpl.ticker.MaxNLocator(prune='both',nbins=4))
						axlist_extent[pnum].set_xlabel('extent '+r'$\mathrm{\sigma_{a,b}\,(nm)}$',
							fontsize=fsaxlabel)
							
					#---modified for IET revisions stage 2014.06.13
					letterlist = ['a','b','c','d','e','f','g','h','i','j','k','l'][3:]
					axlist_extent[pnum].text(0.82,0.86,r'$\mathbf{('+letterlist[pnum].capitalize()+')}$',
						transform=axlist_extent[pnum].transAxes,fontsize=14)							
							
					#---computations
					valid_frames = len(hmaxdat)
					scaled_mean = float(valid_frames)/len(test.data)*mean(hmaxdat)
					result_nums = [
						(round(mean(hmaxdat),4) if valid_frames > 0 else np.nan),
						(round(scaled_mean,4) if valid_frames > 0 else np.nan),
						round(mean(resids),1),
						round(sqrt(mean(sigmas0)),2),round(sqrt(mean(sigmas1)),2),
						round(sqrt(mean(extdat[extdat<2*extent_range])),2),
						round(sqrt(mean(extdat)),2),
						valid_frames]
					print result_nums
					geog = listlook(test.notes,'geography')
					if geog[0] == 'control':
						casual_name = analysis_descriptors[geog[1]]['label_text']
					else: casual_name = nbornames[cnum]
					fp.write(
						str(casual_name).ljust(20)+\
						str((' ' if result_nums[0] > 0 else '')+str(result_nums[0])).ljust(18)+\
						str((' ' if result_nums[1] > 0 else '')+str(result_nums[1])).ljust(18)+\
						str('%1.0f'%valid_frames).ljust(8)+\
						str(str(' ' if result_nums[2] > 0 else '')+str(result_nums[2])).ljust(8)+\
						str(str(' ' if result_nums[3] > 0 else '')+str(result_nums[3])).ljust(10)+\
						str(str(' ' if result_nums[4] > 0 else '')+str(result_nums[4])).ljust(10)+\
						str(str(' ' if result_nums[5] > 0 else '')+str(result_nums[5])).ljust(16)+\
						str(str(' ' if result_nums[6] > 0 else '')+str(result_nums[6])).ljust(25)+'\n')
					#---save the mean Hmax if performing the lattice test for later histogram
					if re.match('lattice.+',dimple_plot_type) and not isnan(result_nums[0]):
						lattice_mean_collect[pnum].append(result_nums[0])
						
					#--save for meta-analysis
					#---only do meta analysis if the params_plot is multidimensional
					if len(shape(params_plot)) > 1:
						dimple_meta_set[anum][gnum] = result_nums

				#---truncate the legend if too large
				ax = axlist[pnum]
				if len(plot_order) < 8: 
					ax.legend(loc='upper left',prop={'size':fsaxlegend_small})
				else:
					which_label = slice(-1,None) if plot_order[0] != 0 else slice(None,1)
					h,l = ax.get_legend_handles_labels()
					ax.legend(h[which_label],l[which_label],loc='upper left',
					prop={'size':fsaxlegend_small})
			
			#---collect the sum of each curve for further analysis
			phists.append(mean(chists,axis=0))
			sysnames = [listlook(params_plot[p][0],'label_text') for p in range(len(params_plot))]

		#---plot settings
		tickwid = 0.5
		for a in range(len(axlist)):
			ax = axlist[a]
			ax.grid(True)
			ax.set_xlim((-hifilt,hifilt))
			ax.set_ylim((0,1.1*maxpeak))
			ax.set_yticklabels([])
			color,proper_name = dimple_colordict_systems(sysnames[a])
			ax.set_ylabel(proper_name,fontsize=fsaxlabel)
			ax.axvline(x=0,ls='--',ymax=1.,ymin=0.,lw=1.5,color='k')
			if a == 0:
				ax.set_title(r'$\textbf{curvature}$',fontsize=16)
			ax.xaxis.set_tick_params(width=tickwid)
			ax.yaxis.set_tick_params(width=tickwid)
		for a in range(len(axlist_extent)):
			ax = axlist_extent[a]
			ax.set_ylim((0,1.1*maxpeak_extent))
			ax.grid(True)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
			ax.set_yticklabels([])
			if a < len(axlist_extent)-1: ax.set_xticklabels([])
			ax.set_xlim((0,extent_range))
			if a == 0:
				ax.set_title(r'$\textbf{extent}$',fontsize=16)
			ax.xaxis.set_tick_params(width=tickwid)
			ax.yaxis.set_tick_params(width=tickwid)

		#--- !!!
		for a in range(len(axlist_areas)):
			ax = axlist_areas[a]
			ax.set_ylim((0,1.1*maxpeak_areas))
			ax.grid(True)
			ax.set_xlim((0,min([len(i.surf) for i in msets])))
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=6))
			ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='upper',nbins=6))
			if a == len(axlist_areas)-1:
				ax.set_xlabel('frame')
			else: 
				ax.set_xlabel('')
				ax.set_xticklabels([])
			ax.set_yticklabels([])
			#---modified for IET revisions stage 2014.06.13
			letterlist = ['a','b','c','d','e','f','g','h','i','j','k','l'][6:]
			ax.text(0.85,0.86,r'$\mathbf{('+letterlist[a].capitalize()+')}$',transform=ax.transAxes,
				fontsize=14)
			ax.xaxis.set_tick_params(width=tickwid)
			ax.yaxis.set_tick_params(width=tickwid)
			if a == len(axlist_areas)-1:
				#---extra legend for areas
				outsidelegendax = ax
				axpos = outsidelegendax.get_position()
				outlegend = outsidelegendax.legend(loc='upper left',prop={'size':12,'weight':'bold'})
				for legobj in outlegend.legendHandles:
					legobj.set_linewidth(3.0)

		for a in range(len(axlist_areas2)):
			ax = axlist_areas2[a]
			#---disabled
			if a == 0 and False:
				ax.set_title(r'$\textbf{areas (+/-)}$')
			ax.grid(True)
			#ax.set_xlim((0,1.1*maxfreq))
			ax.yaxis.tick_right()
			ax.yaxis.set_label_position("right")
			ax.set_ylabel(r'$\textbf{area}\:\textbf{(nm\ensuremath{{}^{2}})}$',
				fontsize=14,rotation=270)
			maxfreq = 0.5
			ax.set_ylim(0,1.1*maxpeak_areas)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=6))
			ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both',nbins=6))
			#---disabled all frequency labels
			if a == len(params_plot)-1:
				ax.set_xlabel('frequency',fontsize=14)
				ax.set_xticklabels([])
			else:
				ax.set_xticklabels([])
			#---modified for IET revisions stage 2014.06.13
			letterlist = ['a','b','c','d','e','f','g','h','i','j','k','l'][9:]
			ax.text(0.75,0.86,r'$\mathbf{('+letterlist[a].capitalize()+')}$',transform=ax.transAxes,
				fontsize=14)
			ax.xaxis.set_tick_params(width=tickwid)
			ax.yaxis.set_tick_params(width=tickwid)

		#---shared area label
		axlist_areas[0].text(1.0,1.03,r'$\textbf{fitted areas}$',ha='center',va='bottom',
			rotation='horizontal',fontsize=16,transform=axlist_areas[0].transAxes)

		#---clean and save
		if dimple_plot_type != 'extrema_meta':
			print '\nstatus: saving plot'
			fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
			plt.savefig(pickles+'fig-dimple3-pubplot-'+bigname+\
				'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
				'filter-'+filter_type+'.'+\
				dimple_plot_type+\
				posneg_select+\
				'.png',dpi=300,bbox_inches='tight')
			if show_plots: plt.show()
			plt.close()
			fp.close()
			print 'status: plot saved'
		
		#---only do meta analysis if the params_plot is multidimensional
		if dimple_plot_type == 'extrema_meta':
			dimple_meta.append(dimple_meta_set)
		
	#---dump meta data
	if dimple_plot_type == 'extrema_meta':
		pickledump(dimple_meta,'fig-dimple3-meta-extrema-'+bigname+'.pkl',directory=pickles)
		
if dimple_plot_type == 'extrema_meta':
	result_col = 1
	result_col_names = ['','-hmax_scaled']
	cutoffs = [100,150,200]
	cutoffs_alpha = [0.35,0.67,1]
	analysis_names_plot = analysis_names[::-1]
	testpairs = [['min_dynamic','max_dynamic'],['min','max']]
	testpairs_titles = ['average','dynamic']
	maxhmax = 0.05
	filter_labels = ['standard','center\nin box','center\nin neighborhood']

	filters = [i['filter_type'] for i in params_plot_settings]
	axeslist = [[[] for j in range(len(testpairs))] for i in range(len(filters))]
	gs = gridspec.GridSpec(len(filters),2,wspace=0.0,hspace=0.0)
	fig = plt.figure()
	for mintest,maxtest in testpairs:
		for filt in filters:
			rown = filters.index(filt)
			coln = testpairs.index([mintest,maxtest])
			ax = fig.add_subplot(gs[rown,coln])
			axeslist[rown][coln] = ax
			for cutoff in cutoffs:
				print cutoff
				#---look up dimple_meta indices for a particular test
				meta_inds_min = [list(array(where([[
					listlook(params_plot[i][j],'geography') == mintest and 
					listlook(params_plot[i][j],'callsign') == aname and
					listlook(params_plot[i][j],'cutoff') == cutoff
					for j in range(len(params_plot[i]))] 
					for i in range(len(analysis_names))])).T[0]) 
					for aname in analysis_names_plot]
				meta_inds_max = [list(array(where([[
					listlook(params_plot[i][j],'geography') == maxtest and 
					listlook(params_plot[i][j],'callsign') == aname and
					listlook(params_plot[i][j],'cutoff') == cutoff
					for j in range(len(params_plot[i]))] 
					for i in range(len(analysis_names))])).T[0]) 
					for aname in analysis_names_plot]
				neg = [dimple_meta[filters.index(filt)][meta_inds_min[i][0]][meta_inds_min[i][1]][result_col] 
					for i in range(len(analysis_names))]
				pos = [dimple_meta[filters.index(filt)][meta_inds_max[i][0]][meta_inds_max[i][1]][result_col] 
					for i in range(len(analysis_names))]	
				ax.plot(range(1,len(neg)+1),pos,'ro-',alpha=cutoffs_alpha[cutoffs.index(cutoff)])
				ax.plot(range(1,len(neg)+1),abs((array(neg))),'bo-',alpha=cutoffs_alpha[cutoffs.index(cutoff)])
			if rown == 0: ax.set_title(testpairs_titles[coln],fontsize=fsaxlabel)
			if coln > 0: ax.set_yticklabels([])
			if rown < len(filters)-1: ax.set_xticklabels([])
			if coln == 0: ax.set_ylabel(r'$\left\langle H_{max}\right\rangle\:\mathrm{}$',fontsize=fsaxlabel)
			if coln == 1: 
				ax.set_ylabel(filter_labels[rown],fontsize=fsaxlabel-4,rotation=270)
				ax.yaxis.set_label_position('right')
			if rown == len(filters)-1:
				ax.set_xticks([i+1 for i in range(len(analysis_names_plot))])
				ax.set_xticklabels([(analysis_descriptors[i])['label'] for i in analysis_names_plot],
					fontsize=fsaxlabel)
				plt.setp(ax.xaxis.get_majorticklabels(),rotation=90)
			ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower',nbins=6))
	for ax in [i for j in axeslist for i in j]:	
		ax.set_xlim((0,len(analysis_names_plot)+1))
		ax.set_ylim((0.,maxhmax))
		ax.grid(True)
	plt.savefig(pickles+'fig-dimple3-meta-extrema-'+bigname+result_col_names[result_col]+'.png',
		dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	
#---additional summary plot, updated and expanded
if 'plot_dimple' in routine and re.match('lattice.+',dimple_plot_type):

	#---compile the plot specifications
	params_plot = compile_plot_params(which_dimple_plot,analysis_descriptors)	

	#---plot the sum of the distributions for clarity (or weight by quantity of valid frames)
	fig = plt.figure(figsize=((1+len(phists))*2,8))
	axeslist = []
	gs = gridspec.GridSpec(len(phists)+2,1,wspace=0.0,hspace=0.0)
	ax = plt.subplot(gs[0:2])
	axeslist.append(ax)
	for p in range(2,len(phists)+2):
		ax = plt.subplot(gs[p])
		axeslist.append(ax)
	
	nnames = [listlook(params_plot[p][0],'label_text') for p in range(len(phists))]
	for p in range(len(phists)):
		color,proper_name = dimple_colordict_systems(nnames[p])
		status('curve = '+str(p+1)+'/'+str(len(params_plot)))
		axeslist[0].plot(hmax_mids,phists[p],lw=2,color=color,
			label=proper_name+' ('+str(round(sum(phists[p]*hmax_mids)/sum(phists[p]),3))+\
			r'$\mathsf{\,(nm^{-1})}$'+')')
		hist,edges = numpy.histogram(lattice_mean_collect[p],
			range=(-hifilt-hist_step/2,hifilt+hist_step/2),bins=hmax_nbins)
		axeslist[p+1].bar(edges[:-1],hist,width=(edges[1]-edges[0]),alpha=0.5,color=color,linewidth=0,
			label=proper_name+' ('+str(round(mean(lattice_mean_collect[p]),3))+r'$\mathsf{\,nm^{-1}}$'+')')
		axeslist[p+1].set_ylabel(proper_name,fontsize=fsaxlabel-2)
	#---plot parameters
	for a in range(len(axeslist)):
		ax = axeslist[a]
		ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')	
		ax.legend(loc='upper left',prop={'size':fsaxlegend_small-(2 if a == 0 else 0)})
		ax.grid(True)
		ax.set_xlim((-hifilt,hifilt))
		ax.set_yticklabels([])
		ax.axvline(x=0,ymax=1.,ymin=0.,lw=1.5,color='k')
		if a == len(axeslist)-1: ax.set_xlabel('$\mathsf{H_{max}\,(nm^{-1})}$',fontsize=16)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)	
		if 0: ax.set_ylabel('distribution of means',fontsize=fsaxlabel)
		if a < len(axeslist)-1: ax.set_xticklabels([])
	axeslist[0].set_ylabel('distribution sum',fontsize=fsaxlabel)
	axeslist[0].set_title('lattice test',fontsize=fsaxtitle)
	
	#---clean and save
	print '\nstatus: saving plot'
	plt.savefig(pickles+'fig-dimple3-'+bigname+\
		'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
		'filter-'+filter_type+'.'+\
		dimple_plot_type+'.summary'\
		'.png',dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	print 'status: plot saved'
	
#---plot summary
if 'plot_radar' in routine:

	#---compile the plot specifications
	params_plot = compile_plot_params('plot_radar',analysis_descriptors)
	
	params_plot = params_plot_master[params_plot_master_names.index('radar')]
	dimple_plot_type = params_plot_master_names[params_plot_master.index(params_plot)]	
	if 'msdats' not in globals():
		#---load the current set of msdats
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---load the pickle to see if it exists	
			filespec = specname_guess(sysname,trajsel)
			pklname = pklprefix+spacetag+filespec+'.pkl'
			msdat = unpickle(pickles+pklname)
			if msdat == None: msdats[anum] = []
			else: msdats[anum] = msdat
	
	#---prepare plot panels	
	fig = plt.figure(figsize=(5,2*len(params_plot)))
	gs = gridspec.GridSpec(len(params_plot),3,wspace=0.0,hspace=0.0)
	gs.update(left=0.0,right=1.0)
	gs2 = gridspec.GridSpec(len(params_plot),1,wspace=0.0,hspace=0.0)
	gs2.update(left=0.75,right=1.0)
	axlist,axlist2,axlist3,axhists = [],[],[],[]
	fp = open(pickles+'calc-dimple3-'+bigname+\
		'.cut'+str(listlook(params_plot[0][0],'cutoff'))+spacetag+\
		'filter-'+filter_type+'.'+\
		dimple_plot_type+'.txt','w')
	
	#---loop over rows
	for pnum in range(len(params_plot)):
		ax = plt.subplot(gs[pnum,0])
		axlist.append(ax)
		ax2 = plt.subplot(gs[pnum,1])
		axlist2.append(ax2)
		ax3 = plt.subplot(gs[pnum,2])
		axlist3.append(ax3)
		
		#---loop over plots within a panel
		for gnum in range(len(params_plot[pnum])):
			expt = params_plot[pnum][gnum]
			aname = listlook(expt,'callsign')
			anum = analysis_names.index(aname)
			mset = msets[anum]
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			nborhoods,nbornames = dimple_generate_neighborhood(expt)
			
			#---lightly plot structure in the background
			ax = axlist[pnum]
			im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,
				lognorm=False,cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=0.65,
				ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.],
				tickskip=int(round(mset.griddims[0]/6,-1)),nine=False)
			
			#---plot dots for the location of the center of each dimple
			dat_ind = dimple_test_lookup(expt,aname)[0]
			maxhxys_inds = msdats[anum][dat_ind].get(['type','maxhxys'])
			target_zones = msdats[anum][dat_ind].get(['type','target_zones'])
			for fr in range(len(maxhxys_inds)):
				abs_pt = target_zones[fr][maxhxys_inds[fr]][:2]
				pt = [abs_pt[i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] for i in range(2)]
				ax.add_patch(plt.Circle((pt[0],pt[1]),
					radius=0.5*mset.griddims[0]/35,color='k',alpha=1.,fc='w',fill=True,lw=1.,ec='k'))
			
			#---plot the centers of the dimple		
			ax2 = axlist2[pnum]
			im = plotter2d(ax2,mset,dat=mean(mset.surf,axis=0)/10.,
				lognorm=False,cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=0.4,
				ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.],
				tickskip=int(round(mset.griddims[0]/6,-1)),nine=False)
			params = msdats[anum][dat_ind].get(['type','params'])
			centers = [[i[2],i[3]] for i in params]
			for fr in range(len(params)):
				abs_pt = centers[fr]
				pt = [abs_pt[i]/mean(mset.vecs,axis=0)[i]*mset.griddims[i] for i in range(2)]
				ax2.add_patch(plt.Circle((pt[0],pt[1]),
					radius=0.5*mset.griddims[0]/35,color='k',alpha=1.,fc='w',fill=True,lw=1.,ec='k'))
			shown_count = sum([all([all([centers[c][i] > 0 and centers[c][j] < mset.vec(c)[j]]) 
				for j in range(2)]) for c in range(len(centers))])
			ax2.text(0.95,0.9,str(shown_count),transform=ax2.transAxes,fontsize=fsaxlabel,
				horizontalalignment='right')

			#---plot the centers of the dimple		
			ax3 = axlist3[pnum]
			im = plotter2d(ax3,mset,dat=mean(mset.surf,axis=0)/10.,
				lognorm=False,cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
				ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
				fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.],
				tickskip=int(round(mset.griddims[0]/6,-1)),nine=False)
			protpts = mean(nborhoods[0],axis=0)
			nnames = [listlook(msdats[anum][0].notes,'neighborhood_name')]
			color = dimple_colordict(nnames[0],listing=nnames,scheme=dimple_plot_type,nprots=nprots)
			plothull(ax3,protpts,mset=mset,subdivide=None,
				c=color,alpha=1.,fill=False)
				
			#---plot settings
			ax3.set_ylabel(listlook(expt,'label'),rotation=270,fontsize=fsaxlabel)
			ax3.yaxis.set_label_position("right")
		
		#---plot settings
		for a in range(len(axlist)):
			if a < len(axlist)-1: 
				axlist[a].set_xticklabels([])
				axlist[a].set_xlabel('')
			if a == 0: axlist[a].set_title(r'$H_{max}(x,y)$',fontsize=fsaxlabel)
		for a in range(len(axlist2)):
			axlist2[a].set_yticklabels([])
			axlist2[a].set_ylabel('')
			if a < len(axlist)-1: 
				axlist2[a].set_xticklabels([])
				axlist2[a].set_xlabel('')
			if a == 0: axlist2[a].set_title(r'dimple centers',fontsize=fsaxlabel)
		for a in range(len(axlist3)):
			axlist3[a].set_yticklabels([])
			if a < len(axlist)-1: 
				axlist3[a].set_xticklabels([])
				axlist3[a].set_xlabel('')
				
	#---clean and save
	print '\nstatus: saving plot'
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	plt.savefig(pickles+'fig-dimple3-'+bigname+\
		'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
		'filter-'+filter_type+'.'+\
		dimple_plot_type+\
		'.png',dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	fp.close()
	print 'status: plot saved'
	
#---plot summary
if 'plot_topog' in routine:
	#---compile the plot specifications
	params_plot = compile_plot_params(topog_plot_type,analysis_descriptors)	
	if 'msdats' not in globals():
		#---load the current set of msdats
		msdats = [[] for i in range(len(analysis_names))]
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			anum = analysis_names.index(aname)
			#---load the pickle to see if it exists	
			filespec = specname_guess(sysname,trajsel)
			pklname = pklprefix+spacetag+filespec+'.pkl'
			msdat = unpickle(pickles+pklname)
			if msdat == None: msdats[anum] = []
			else: msdats[anum] = msdat
	topog_nbins = 30
	params_plot = compile_plot_params(topog_plot_type,analysis_descriptors)
	dats = [[array(flatten(msdats[i][k].data)) for k in range(len(msdats[i]))] for i in range(len(msdats))]
	maxz = round(max([abs(min(flatten(dats))),abs(max(flatten(dats)))]),-1)+10
	extremz = max([max([mean(mset.surf,axis=0).max(),abs(mean(mset.surf,axis=0).min())])
		for mset in msets])
	fig = plt.figure(figsize=((1+len(dats))*3.75,12))
	axeslist,structlist = [],[]
	gssnap = gridspec.GridSpec(5,2*len(dats),wspace=0.2,hspace=0.2)
	gssnap.update(left=0.0,right=1.0,top=1.0,bottom=0.7)
	gs = gridspec.GridSpec(1,len(dats),wspace=0.0,hspace=0.0)
	gs.update(left=0.0,right=1.0,top=0.65,bottom=0.35)
	gs2 = gridspec.GridSpec(1,len(dats),wspace=0.0,hspace=0.0)
	gs2.update(left=0.0,right=1.0,top=0.3,bottom=0.0)
	nnames = [listlook(msdats[0][i].notes,'neighborhood_name') for i in range(len(params_plot[0]))]
	maxfreq = 0
	for pnum in range(len(params_plot)):
		ax = plt.subplot(gs2[pnum])
		axeslist.append(ax)
		for gnum in range(len(params_plot[pnum])):
			expt = params_plot[pnum][gnum]
			aname = listlook(expt,'callsign')
			anum = analysis_names.index(aname)
			mset = msets[anum]
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			nborhoods,nbornames = dimple_generate_neighborhood(expt)
			dat_ind = dimple_test_lookup(expt,aname)[0]
			hist,edges = numpy.histogram(flatten(msdats[anum][dat_ind].data),range=(-maxz,maxz),
				bins=topog_nbins)
			mids = 1./2*(edges[1:]+edges[:-1])/10.
			color = dimple_colordict(nnames[gnum],listing=nnames,scheme=topog_plot_type)
			ls = '--' if nnames[gnum] in ['valley (dynamic)','peak (dynamic)'] else '-'
			axeslist[pnum].plot(mids,hist,ls,c=color,
				label=listlook(msdats[anum][dat_ind].notes,'neighborhood_name'))
			maxfreq = max([maxfreq,max(hist)])
		ax.legend(loc='upper left',fontsize=fsaxlegend-4)
		ax.set_xlabel(r'$z_{max},z_{min}\:\mathrm{(nm)}$',fontsize=fsaxlabel)
		ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
		if pnum == 0: ax.set_ylabel('extrema (dynamic)',fontsize=fsaxlabel)
		ax = plt.subplot(gs[pnum])
		structlist.append(ax)
		im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0)/10.,
			lognorm=False,cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=0.65,
			ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
			fs=fsaxlabel,label_style='xy',lims=[-extremz/10.,extremz/10.],
			tickskip=int(round(mset.griddims[0]/6,-1)),nine=False)
		if nprots > 0:
			protpts = mean(mset.protein,axis=0)
			plothull(ax,protpts,mset=mset,subdivide=nprots,
				c='k',alpha=1.,fill=False)
		if pnum == len(params_plot)-1:
			axins2 = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins2,orientation="vertical")
			axins2.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{nm})$',
				fontsize=fsaxlabel,rotation=270)
			axins2.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='both'))
		#---top row of birdseye images
		img = mpimg.imread(pickles+imagelist[aname])
		ax = plt.subplot(gssnap[2:5,2*pnum])
		imgplot = ax.imshow(img)
		ax.set_yticks([])
		ax.set_xticks([])
		shifter = [100,-100] if aname == 'v614-120000-220000-200' else [0,0]
		ax.add_patch(plt.Rectangle((800-shifter[0],800-shifter[1]),
			width=1400,height=1400,fill=False,ec='y',lw=2))
		#---top row of birdseye images
		img = mpimg.imread(pickles+imagelist[aname])
		ax = plt.subplot(gssnap[2:5,2*pnum+1])
		imgplot = ax.imshow(img[800+shifter[0]:2200+shifter[0],800+shifter[1]:2200+shifter[1]])
		ax.set_yticks([])
		ax.set_xticks([])
		#---second row of birdseye images
		img = mpimg.imread(pickles+imagelist2[aname])
		ax = plt.subplot(gssnap[0:2,2*pnum:2*pnum+2])
		ax.set_title(label,fontsize=fsaxtitle)
		shifter = -150 if aname == 'v612-75000-175000-200' else 0
		imgplot = ax.imshow(img[500+shifter:1500+shifter,700:3300])
		ax.set_yticks([])
		ax.set_xticks([])
	for ax in axeslist:
		ax.set_ylim((0,1.1*maxfreq))
		ax.grid(True)
		ax.set_yticklabels([])
		ax.set(aspect=float(maxz*2/maxfreq/10))
	for a in range(len(structlist)):
		ax = structlist[a]
		if a > 0 and 0: ax.set_ylabel('')
		if a > 0 and 0: ax.set_yticklabels([])
	
	#---clean and save
	print '\nstatus: saving plot'
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	plt.savefig(pickles+'fig-topography-summary-'+bigname+\
		'.cut'+str(listlook(expt,'cutoff'))+spacetag+\
		('filter-'+filter_type+'.' if 'filter_type' in globals() else '')+\
		dimple_plot_type+\
		'.png',dpi=300,bbox_inches='tight')
	if show_plots: plt.show()
	plt.close()
	print 'status: plot saved'	
	

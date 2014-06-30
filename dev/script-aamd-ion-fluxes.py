#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = 'name P'

analysis_descriptors = {
	'v514-10000-29000-100':
		{'sysname':'membrane-v514',
		'sysname_lookup':'membrane-v514-ions',
		'trajsel':'s3-sim-compbio-md.part0004.10000-29000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v514.a2-surfacer.s3-sim-compbio-md.part0004.10000-29000-100.pkl',
		'ionname':'NA'},
	'v532-20000-58000-100':
		{'sysname':'membrane-v532',
		'sysname_lookup':'membrane-v532-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-58000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',
		'ionname':'Cal'},
	'v531-20000-62000-100':
		{'sysname':'membrane-v531',
		'sysname_lookup':'membrane-v531-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-62000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v531.a6-surfacer.s4-sim-trestles-md.part0007.20000-62000-100.pkl',
		'ionname':'MG',
		'sysname_lookup_struct':'membrane-v531-atomP',
		'trajsel_struct':'s4-sim-trestles-md.part0007.20000-62000-100.atomP.xtc',
		'director':director_aamd_asymmetric},
	'v530-30000-100000-100':
		{'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-ions',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl',
		'ionname':'NA'},
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl',
		'ionname':'Cal'},
	'v509-40000-90000-10':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-binding',
		'trajsel':'s4-kraken-md.part0007.20000-90000-10.binding.xtc',
		'sysname_lookup_struct':'membrane-v509-binding',
		'trajsel_struct':'s4-kraken-md.part0007.20000-90000-10.binding.xtc',
		'ionname':'NA',
		'director':director_aamd_symmetric}}
analysis_names = [
	'v532-20000-58000-100',
	'v530-30000-100000-100',
	'v531-20000-62000-100',
	'v509-40000-90000-10'
	][-1:]
routine = ['load','compute_z','compute_radial_binary'][-1:]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def bilayer_shift(ionspos=None,midz=None,mset=None,vecs=None):
	'''Function to change coordinate systems so the midplane is at zero.'''
	#---ionspos holds absolute ion positions
	#---midz holds the average midplane z-position for each frame
	if ionspos == None: ionspos = globals()['ionspos']
	if mset == None: mset = globals()['mset']
	if midz == None: midz = globals()['midz']
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	deltaz = array(ionspos)[...,2]-tile(midz,(shape(ionspos)[1],1)).T
	abovez = (deltaz>tile(vecs[:,2]/2.,(shape(deltaz)[1],1)).T)
	belowz = (deltaz<-1*tile(vecs[:,2]/2.,(shape(deltaz)[1],1)).T)
	relpos = deltaz - abovez*tile(vecs[:,2],(shape(abovez)[1],1)).T + \
		belowz*tile(vecs[:,2],(shape(belowz)[1],1)).T
	return relpos

def binlister(method,monoz=None,mset=None,binw=5,custom_posts=None):
	'''This function generates a list of bin limits for a particular system. Feed this to discretizer_z.'''
	if mset == None: mset = globals()['mset']
	if monoz == None: monoz = globals()['monoz']	
	if method == 'outside-inside':
		#---this method uses a set of bins inside the bilayer and outside, so the monolayer definitions are 
		#---...hard-coded as bin edges. bins are determined frame-wise
		#---midplane positions
		midz = mean(monoz,axis=1)
		#---box vectors
		vecs = [mset.vec(i) for i in range(mset.nframes)]
		#---average distance outside of phosphate groups
		waterdist = mean(vecs,axis=0)[2]-mean([i[0]-i[1] for i in monoz])
		#---bin counts in the water
		nbins_out = int(round(waterdist/binw))
		nbins_in = int(round(mean([i[0]-i[1] for i in monoz])/binw))
		#---generate list of bin zones
		binlists = []
		for fr in range(mset.nframes):
			thick_in = (monoz[fr,0]-monoz[fr,1])/2.
			thick_out = (vecs[fr][2]-thick_in*2)/2.
			if (nbins_out % 2) == 1:
				end_binw = thick_out/(nbins_out/2+0.5)
				thick_out = thick_out - end_binw
				binlist = sum([list(i) for i in [
					linspace(-thick_out-end_binw-thick_in,-thick_out,1+1)[:-1],
					linspace(-thick_in-thick_out,-thick_in,nbins_out/2+1)[:-1],
					linspace(-thick_in,thick_in,nbins_in+1)[:-1],
					linspace(thick_in,thick_in+thick_out,nbins_out/2+1)[:-1],
					linspace(thick_in+thick_out,thick_in+thick_out+end_binw,2+1)[:]]])
			else: 
				binlist = sum([list(i) for i in [
					linspace(-thick_in-thick_out,-thick_in,nbins_out/2+1)[:-1],
					linspace(-thick_in,thick_in,nbins_in+1)[:-1],
					linspace(thick_in,thick_in+thick_out,nbins_out/2+1)[:]]])
			binlists.append(binlist)
	elif method == 'fixed':
		#---midplane positions
		midz = mean(monoz,axis=1)
		#---box vectors
		vecs = [mset.vec(i) for i in range(mset.nframes)]
		#---generate bin edges based on means first in order to figure out bin counts
		relvec = mean(vecs,axis=0)[2]/2.
		relmonoz = mean([monoz[fr]-midz[fr] for fr in range(mset.nframes)],axis=0)
		meanbinlist = [
			arange(relmonoz[1],-relvec-binw,-binw)[:0:-1],
			arange(relmonoz[1],0,binw),
			arange(relmonoz[0],0,-binw)[::-1],
			arange(relmonoz[0],relvec+binw,binw)[1:]]
		#---generate list of bin zones
		binlists = []
		for fr in range(mset.nframes):
			relvec = vecs[fr][2]/2.
			relmonoz = monoz[fr]-midz[fr]
			#---generate the bins for a simulation segment, then consolidate to match the counts of the mean
			bl1 = arange(relmonoz[1],-relvec-binw,-binw)[:0:-1]
			bl2 = arange(relmonoz[1],0,binw)
			bl3 = arange(relmonoz[0],0,-binw)[::-1]
			bl4 = arange(relmonoz[0],relvec+binw,binw)[1:]
			if len(bl1) > len(meanbinlist[0]): bl1 = concatenate((bl1[0:1],bl1[2:]))
			elif len(bl1) < len(meanbinlist[0]): bl1 = concatenate(([bl1[0]-binw],bl1))
			if len(bl2) > len(meanbinlist[1]): bl2 = concatenate((bl2[0:-2],bl2[-1:]))
			elif len(bl2) < len(meanbinlist[1]): bl2 = concatenate((bl2,[0.]))
			if len(bl3) > len(meanbinlist[2]): bl3 = concatenate((bl3[0:1],bl3[3:]))
			elif len(bl3) < len(meanbinlist[2]): bl3 = concatenate(([0.],bl3))
			if len(bl4) > len(meanbinlist[3]): bl4 = concatenate((bl4[0:-2],bl4[-1:]))
			elif len(bl4) < len(meanbinlist[3]): bl4 = concatenate((bl4,[bl4[-1]+binw]))
			binlist = concatenate((bl1,bl2,bl3,bl4))
			binlists.append(list(binlist))
	elif method == 'custom':
		#---midplane positions
		midz = mean(monoz,axis=1)
		#---box vectors
		vecs = [mset.vec(i) for i in range(mset.nframes)]
		binlists = []
		for fr in range(mset.nframes):
			relvec = vecs[fr][2]/2.
			relmonoz = monoz[fr]-midz[fr]
			binlist = \
				[-relvec-binw]+[relmonoz[1]-i for i in custom_posts][::-1] +\
				[relmonoz[0]+i for i in custom_posts]+[relvec+binw]
			binlists.append(binlist)
	else: 
		raise Exception('except: you must specify a binning scheme')
	#---calculate bins corresponding to monolayers
	monobins = [[where((monoz[i,j]-midz[i])<binlists[i][1:])[0][0]+[0,-1][j] 
		for i in range(len(monoz))] for j in range(2)]
	return binlists,[int(i) for i in list(mean(monobins,axis=1))]

#---discretize the trajectory

def discretizer_z(binedges,relpos=None):
	'''Discretizes an ion trajectory in the z-dimension given a set of positions and a list of bins.'''
	if relpos == None: relpos = globals()['relpos']
	disctraj  = array([[where(relpos[j,i]<binedges[j][1:])[0][0] for i in range(len(relpos[j]))] 
		for j in range(len(binedges))]).T
	return disctraj
	
def discretizer_radial_binary(selstring,ionname,mset=None,frameslice=None):
	'''Discretizes a trajectory into inside-outside bins for a radial distance defintion.'''
	if mset == None: mset = globals()['mset']
	if frameslice == None: frameslice = slice(None,None)
	allresids = sort(mset.universe.selectAtoms('name '+ionname).resids())
	ts = range(mset.nframes)[frameslice]
	starttime = time.time()
	disctrajt = []
	for t in ts:
		status('t = '+str(t),start=starttime,i=t,looplen=len(ts))
		mset.gotoframe(t)
		sel = mset.universe.selectAtoms(selstring)
		inbin = [where(allresids==i)[0][0] for i in sel.resids()]
		discframe = zeros(len(allresids))
		discframe[inbin] = 1
		disctrajt.append(discframe)
	disctraj = array(disctrajt).T	
	return disctraj

#---define zones for plotting calculations relative to discrete trajectories

def define_zones(zonetype):
	if zonetype == 'radial_binary':
		zonelist = [0,1]
		zonecolors = ['b','r']
		zonealphas = [0.5,0.5]
		zonelabels = ['free','bound']
		return [zonelist,zonecolors,zonealphas,zonelabels]
		
#---plotting and analysis functions for the discrete trajectories

def plot_ion_distribution(disctraj,binedges,fig=None,thisax=None,text=True,colors=True,mset=None):
	'''Generic plotting function for ion distributions by height from midplane according to binedges.'''
	#---Nb only works with discretizer_z
	if fig == None: fig = plt.figure()
	if mset == None: mset = globals()['mset']
	gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
	if thisax == None: ax = fig.add_subplot(gs[0])
	else: ax = thisax
	counts,edges = histogram(disctraj.flatten(),range=(0,len(binedges[0])-1),
		bins=len(binedges[0])-1,normed=True)
	meanedges = mean(binedges,axis=0)
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	lefts = meanedges[:-1]	
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	halfvec = mean(vecs,axis=0)[2]/2
	for n in range(len(mids)):
		if mids[n] == 0: color = 'k'
		elif mids[n] > 0: color = 'r'
		elif mids[n] < 0: color = 'b'
		ax.bar(lefts[n],counts[n],color=color,width=list(meanedges[1:]-meanedges[:-1])[n],
			alpha=1-abs(mids[n]/halfvec))
	ax.grid(True)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_xlabel(r'$z(x,y)\,\mathrm{(nm)}$',fontsize=fsaxlabel)
	ax.set_ylabel('frequency',fontsize=fsaxlabel)
	if text == False:
		ax.set_ylabel('')
		ax.set_xlabel('')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
	plt.show()
	return counts,edges
	
def plot_ion_transition_matrix(disctraj,binedges,fig=None):
	'''Plot the transition matrix for the discrete ion trajectories.'''
	#---Nb only works with discretizer_z
	if fig == None: fig = plt.figure()
	gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
	ax = fig.add_subplot(gs[0])
	#---prepare the transition matrix
	disctrajt = disctraj.T
	tmat = zeros((len(binedges[0])-1,len(binedges[0])-1))
	for fr in range(len(disctrajt)-1):
		for i in range(len(disctrajt[fr])):
			tmat[disctrajt[fr][i],disctrajt[fr+1][i]] += 1
	#---plot the transition matrix
	plt.imshow(tmat,interpolation='nearest',origin='lower',
		cmap=mpl.cm.binary)
	tickstep = 2
	meanedges = mean(binedges,axis=0)
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	#---indexing to number about zero in steps
	tickspots = [i+len(mids)/2 for i in range(0,-len(mids)/2,-tickstep)[:0:-1]+range(0,len(mids)/2,tickstep)]
	ax.set_xticks(tickspots)
	ax.set_yticks(tickspots)
	ax.set_yticklabels([int(round(mids[i])) for i in tickspots])
	ax.set_xticklabels([int(round(mids[i])) for i in tickspots])
	ax.set_title('ion transition matrix',fontsize=fsaxtitle)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_xlabel(r'$z_t(x,y)\,\textrm{\AA}$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$z_{t+1}(x,y)\,\textrm{\AA}$',fontsize=fsaxlabel)
	plt.show()
	
def plot_ion_residence_z(disctraj,binedges,fig=None,ignorezero=True,maxtime=100,
	mset=None,dofit=False,cumulative=True):
	#---Nb only works with discretizer_z
	'''Plot residence time distributions.'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	#---fineness of the time bins is currently hard-coded
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0],[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	#---this step implicitly filters out ions that stay bound the entire time
	maxresid = max([max(i) for i in jumpdurations if len(i) != 0])
	jumphists = [histogram(jumpdurations[i],range=(0,ntimes),bins=ntimes+1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	residences = [[] for i in range(len(binedges[0])-1)]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[jumpbins[k][i]].append(jumpdurations[k][i])
	maxresid = max([max(i) for i in residences if i != []])
	times = mset.time_dt*histogram(residences[0],range=(0,maxresid),bins=ntimes)[1][1:]
	resdists = [histogram(residences[i],range=(0,maxresid),bins=ntimes)[0] 
		for i in range(len(residences)) if i != []]	
	#---plot
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	halfvec = mean(vecs,axis=0)[2]/2
	mean(vecs,axis=0)[2]
	meanedges = mean(binedges,axis=0)
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	ax = plt.subplot(111)
	for r in range(len(resdists)):
		if ignorezero == False or mids[r] != 0:
			dat = resdists[r]
			indscut = where((1*(times<maxtime)+1*(dat!=0))==2)[0]
			if mids[r] > 0.: c = 'r'
			elif mids[r] == 0: c = 'k'
			else: c = 'b'
			inds = where(dat!=0)[0]
			if len(dat[indscut]) > 0 and dofit:
				[tau,intercept] = numpy.polyfit(log(times[indscut]),log(dat[indscut]),1)
				fitted = array([exp(intercept)*exp(tau*log(i)) for i in times])
				ax.plot(times[inds],fitted[inds],'--',lw=0.5,c=c)
			else: tau = 0
			if cumulative:
				ax.plot(times[inds],sum(dat[inds])-cumsum(dat[inds]),'o-',c=c,mec=c,lw=1,
					alpha=1-abs(mids[r]/halfvec),
					label=str(mids[r]/10.)+'('+str(round(-tau,2))+')')
			else:
				ax.plot(times[inds],dat[inds],'o-',c=c,mec=c,lw=1,alpha=1-abs(mids[r]/halfvec),
					label=str(mids[r]/10.)+'('+str(round(-tau,2))+')')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((1,ntimes*mset.time_dt*2))
	ax.set_xlabel(r'$\tau\,(ps)$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$frequency$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_title('residence times',fontsize=fsaxtitle)
	ax.legend(loc='upper right')
	ax.grid(True)
	plt.show()

def plot_ion_time_correlations(disctraj,binedges,fig=None,
	thisax=None,mset=None,t0s_select=None,showlegend=False,plotall=False):
	'''Plot the time correlation function for presence in each zone.'''
	#---Nb only works with discretizer_z
	if fig == None: fig = plt.figure(figsize=(5,6))
	if mset == None: mset = globals()['mset']
	disctrajt = disctraj.T
	meanedges = mean(binedges,axis=0)
	mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
	times = mset.time_dt*arange(mset.nframes)+10.
	nslices = len(binedges[0])-1
	t0s = range(0,len(disctrajt)-1)
	starttime = time.time()
	init_occupants = []
	for t0 in t0s:
		init_occupants_by_bin = []
		for binn in range(nslices):
			init_occupants_by_bin.append(where(disctrajt[t0]==binn)[0])
		status('t0 = '+str(t0),start=starttime,i=t0,looplen=len(t0s))
		init_occupants.append(init_occupants_by_bin)
	if thisax == None: ax = plt.subplot(111)
	else: ax = thisax
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	halfvec = mean(vecs,axis=0)[2]/2
	starttime = time.time()
	if t0s_select == None:
		t0s = range(0,len(disctrajt)-1)
	else: t0s = t0s_select
	collect_rts = [[] for i in range(nslices)]
	for t0 in t0s:
		rt = [[(float(sum(disctrajt[i][init_occupants[t0][binn]]==binn))/len(init_occupants[t0][binn]) 
			if len(init_occupants[t0][binn]) > 0 else []) for i in range(t0,len(disctrajt))] 
			for binn in range(nslices)]
		if mean([i[0] for i in rt if i[0] != []]) != 1.:
			raw_input('status: calculation failure')
		for n in range(nslices):
			if mids[n] > 0.: c = 'r'
			elif mids[n] == 0: c = 'k'
			else: c = 'b'
			alpha = (1.-(mids[n]-2)/5.)
			if mids[n] != 0. and all([type(i)==float for i in rt[n]]):
				if plotall == True:
					ax.plot(times[:len(rt[n])],rt[n],c=c,alpha=alpha)
				collect_rts[n].extend(list(array([times[:len(rt[n])],rt[n]]).T))
		status('t0 = '+str(t0),start=starttime,i=t0,looplen=len(t0s))
	collect_rts = [array(i) for i in collect_rts]
	for n in range(nslices):
		if collect_rts[n] != []:
			tmp = collect_rts[n]
			if mids[n] > 0.: c = 'r'
			elif mids[n] == 0: c = 'k'
			else: c = 'b'
			timelist = list(sort(list(set(tmp[:,0]))))
			avgcurv = [mean(tmp[:,1][tmp[:,0]==i]) for i in timelist]
			ax.plot(timelist,avgcurv,lw=2,c=c,
				alpha=1-abs(mids[n]/halfvec),
				label='('+str(round(meanedges[n]/10.,1))+','+str(round(meanedges[n+1]/10.,1))+')')
	#---plot settings
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim((10**-1,2))
	ax.set_xlim((10**0,max(times)*10))
	ax.grid(True)
	ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$\mathrm{R(\mathbf{\tau})\:(ps)}$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_title('time correlations',fontsize=fsaxtitle)
	if showlegend: ax.legend(loc='lower right')
	#---inset ion distribution
	axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="30%",height="30%",loc=3)
	plot_ion_distribution(disctraj,binedges,thisax=axins,fig=fig,text=False,colors=True,mset=mset)
	plt.savefig(pickles+'fig-ion_time_correlate-'+specname_guess(sysname,trajsel).strip('membrane-')+'.png',
		dpi=500,bbox_inches='tight')
	plt.show()
	
def plot_ion_temporary_experiment():
	'''Testing code.'''
	if 1:
		fig=None
		ignorezero=True
		maxtime=100;
		mset = mset_surf
		dofit=False
		cumulative=False
		'''Plot residence time distributions.'''
		spec = 3
		if fig == None: fig = plt.figure(figsize=(8,10))
		if mset == None: mset = globals()['mset']
		ntimes = len(disctraj[0])
		jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0],[len(disctraj[0])])) 
			for i in range(len(disctraj))]
		jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
		jumphists = [histogram(jumpdurations[i],range=(1,ntimes),bins=ntimes-1)[0] 
			for i in range(len(jumpdurations))]
		jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
		residences = [[] for i in range(len(binedges[0])-1)]
		for k in range(len(jumpbins)):
			if len(jumpbins[k]) > 0:
				for i in range(len(jumpbins[k])):
					residences[jumpbins[k][i]].append(jumpdurations[k][i])
		times = mset.time_dt*histogram(residences[0],range=(1,ntimes),bins=ntimes-1)[1][:-1]
		resdists = [histogram(residences[i],range=(1,ntimes),bins=ntimes-1)[0] 
			for i in range(len(residences)) if i != []]	

		#---plot
		vecs = array([mset.vec(i) for i in range(mset.nframes)])
		halfvec = mean(vecs,axis=0)[2]/2
		mean(vecs,axis=0)[2]
		meanedges = mean(binedges,axis=0)
		mids = [int(round(i)) for i in 1./2*(meanedges[1:]+meanedges[:-1])]
		ax = plt.subplot(111)
		for r in range(len(resdists)):
			if ignorezero == False or mids[r] != 0:
				dat = resdists[r]
				if mids[r] > 0.: c = 'r'
				elif mids[r] == 0: c = 'k'
				else: c = 'b'
				ax.plot(times,dat,'o-',c=c,mec=c,lw=1,alpha=1-abs(mids[r]/halfvec),
					label=str(mids[r]/10.))
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlim((1,ntimes*mset.time_dt*2))
		ax.set_xlabel(r'$\tau\,(ps)$',fontsize=fsaxlabel)
		ax.set_ylabel(r'$frequency$',fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		ax.set_title('residence times',fontsize=fsaxtitle)
		ax.legend(loc='upper right')
		ax.grid(True)
		plt.show()
	if 0:
		ax = plt.subplot(111)
		valid_inds = where(array([len(i) for i in alldat])>0)[0]
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlim((10**0,max(times)*10))
		ax.plot(times[valid_inds],[sum(alldat[i]) for i in valid_inds],'o-')
		plt.show()
	if 1:
		bin_ind = 3
		ax = plt.subplot(111)
		alldat = [[] for i in range(4999)]
		for ion_ind in range(len(jumpdurations)):
			counts = histogram(jumpdurations[ion_ind][jumpbins[ion_ind]==bin_ind],
				range=(1,ntimes),bins=ntimes-1)[0]
			for k in where(counts!=0)[0]:
				alldat[k].append(counts[k])
			ax.plot(times[counts!=0],times[counts!=0]*counts[counts!=0]/\
				sum(times[counts!=0]*counts[counts!=0]),'-')
		ax.plot(times,times*mean(alldat,axis=0)/sum(times*mean(alldat,axis=0)),'o-',lw=2,c='k')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlim((10**0,max(times)*10))
		ax.grid(True)
		ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
		ax.set_ylabel(r'$\mathrm{R(\mathbf{\tau})\:(ps)}$',fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		ax.set_title('time correlations',fontsize=fsaxtitle)
		plt.show()
		
def plot_ion_residence_radial_binary(disctraj,zonetype,fig=None,ignorezero=True,
	mset=None,dofit=False,cumulative=False,scale_by_time=None):
	'''Plot residence time distributions, radial binary version.'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0]+1,[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	jumphists = [histogram(jumpdurations[i],range=(0,ntimes),bins=ntimes+1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	#---get all zone definitions for plotting
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	residences = [[] for i in range(len(zonelist))]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[int(jumpbins[k][i])].append(jumpdurations[k][i])
	times = mset.time_dt*histogram(residences[0],range=(1,ntimes+1),bins=ntimes)[1][:-1]
	resdists = [histogram(residences[i],range=(1,ntimes+1),bins=ntimes)[0] 
		for i in range(len(residences)) if i != []]
	#---plot
	ax = plt.subplot(111)
	extrem = [10**30,0]
	for r in range(len(resdists)):
		dat = resdists[r]
		c = zonecolors[r]
		alpha = zonealphas[r]
		label = zonelabels[r]
		inds = where(dat!=0)[0]
		if scale_by_time == True:
			if cumulative == False:
				plotdat = dat[inds]*times[inds]/float(times[-1])/len(disctraj)
			else:
				plotdat = cumsum(dat[inds]*times[inds]/float(times[-1])/len(disctraj))
			ax.plot(times[inds],plotdat,
				'o-',c=c,mec=c,lw=1,alpha=alpha,label=label+' ('+\
				str(round(float(sum(disctraj==zonelist[r]))/product(shape(disctraj)),2))+
				(','+str(int(plotdat[-1]*len(disctraj)))+'bound' 
					if (times[inds][-1] == mset.time_dt*ntimes and cumulative == False) else '')+')')
		else:
			plotdat = dat[inds]
			ax.plot(times[inds],plotdat,'o-',c=c,mec=c,lw=1,alpha=alpha,label=label)
		if extrem[0] > min(plotdat): extrem[0] = min(plotdat)
		if extrem[1] < max(plotdat): extrem[1] = max(plotdat)
	if scale_by_time == True: ax.set_ylim((0.5*extrem[0],10*extrem[1]))
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim((1,ntimes*mset.time_dt*2))
	ax.set_xlabel(r'$\tau\,(ps)$',fontsize=fsaxlabel)
	if scale_by_time == True: ax.set_ylabel(r'$\mathrm{P(\tau)}$',fontsize=fsaxlabel)
	else: ax.set_ylabel(r'$frequency$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.legend(loc='upper right')
	ax.grid(True)
	if cumulative and scale_by_time: extratag = '-cdf'
	elif scale_by_time and not cumulative: extratag = '-pdf'
	else: extratag = ''
	plt.savefig(pickles+'fig-ion_residence-'+specname_guess(sysname,trajsel).strip('membrane-')+\
		extratag+'.png',dpi=500,bbox_inches='tight')
	plt.show()
	
def plot_ion_residence_radial_binary_dev(disctraj,zonetype,fig=None,
	mset=None,normed=True,scale_by_time=None,residence_exact=False):
	'''Plot residence time distributions, radial binary version.'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	#---header is the same as the original residence time calculator
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0]+1,[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	jumphists = [histogram(jumpdurations[i],range=(0,ntimes),bins=ntimes+1)[0] 
		for i in range(len(jumpdurations))]
	jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	residences = [[] for i in range(len(zonelist))]
	for k in range(len(jumpbins)):
		if len(jumpbins[k]) > 0:
			for i in range(len(jumpbins[k])):
				residences[int(jumpbins[k][i])].append(jumpdurations[k][i])
	times = mset.time_dt*histogram(residences[0],range=(1,ntimes+1),bins=ntimes)[1][:-1]
	resdists = [histogram(residences[i],range=(1,ntimes+1),bins=ntimes)[0] 
		for i in range(len(residences)) if i != []]
	zonetype = 'radial_binary'
	times = mset.time_dt*arange(mset.nframes)+10.
	zonelist,zonecolors,zonealphas,zonelabels = define_zones(zonetype)
	disctrajt = array(disctraj).T
	t0s = range(0,len(disctrajt)-1)
	starttime = time.time()
	init_occupants = []
	for t0 in t0s:
		init_occupants_by_bin = []
		for binn in zonelist:
			init_occupants_by_bin.append(where(disctrajt[t0]==binn)[0])
		init_occupants.append(init_occupants_by_bin)
	collected = [[] for i in zonelist]
	collected_residence = [[] for i in zonelist]
	for t0 in t0s:
		for z in zonelist:
			inds = init_occupants[t0][z]
			staytimes = array([jumptimes[i][jumptimes[i]>t0][0]-t0 for i in inds])
			counts,edges = histogram(staytimes,range=(1,ntimes+1),bins=ntimes)
			collected[z].append(counts)
			if residence_exact:
				residence_curve = [sum(staytimes>i) for i in range(ntimes-t0)]#+[0 for i in range(t0)]
				collected_residence[z].append(residence_curve)
		status('t0 = '+str(t0),start=starttime,i=t0,looplen=len(t0s))
	if residence_exact:
		meanrescurves = [[mean([collected_residence[k][t0][j] 
			for t0 in t0s if j < len(collected_residence[k][t0])]) 
			for j in range(ntimes)] 
			for k in range(len(collected_residence))]
	ax = plt.subplot(111)
	for z in zonelist:
		if normed and not residence_exact:
			#---could normalize by len(mean(collected[z],axis=0)) or sum(mean(collected[z],axis=0)) also
			ax.plot(times,mean(collected[z],axis=0)/mean(collected[z],axis=0)[0],lw=2,
				c=zonecolors[z],label=zonelabels[z])
		elif not normed and not residence_exact:
			ax.plot(times,mean(collected[z],axis=0),lw=2,
				c=zonecolors[z],label=zonelabels[z])
		elif residence_exact:
			ax.plot(times,meanrescurves[z],'.-',lw=1,
				c=zonecolors[z],label=zonelabels[z])
	ax.set_xlim((1,2*times[-1]))
	ax.set_xscale('log')
	#ax.set_yscale('log')
	ax.grid(True)
	ax.set_ylabel(r'$\mathrm{N_0}$',fontsize=fsaxlabel)
	ax.set_xlabel(r'$\mathrm{\mathbf{\tau}\:(ps)}$',fontsize=fsaxlabel)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.legend(loc='lower left')
	if residence_exact: extratag = '_dev_exact'
	elif normed: extratag = '_dev_normed'
	else: extratag = '_dev'
	extratag2 = '-'+str(buffer_size)
	plt.savefig(pickles+'fig-ion_residence'+extratag+'-'+\
		specname_guess(sysname,trajsel).strip('membrane-')+extratag2+'.png',
		dpi=500,bbox_inches='tight')
	plt.show()

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---generic load routine for analyzing ions relative to the bilayer position
if 'load' in routine or 'ionspos' not in globals():
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		trajsel = trajsel_struct
		sysname_lookup = sysname_lookup_struct
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			keysysname='sysname_lookup_struct',keytrajsel='trajsel_struct')
		mset_surf = MembraneSet()
		mset_surf.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		#---collect ion positions
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		#---collect monolayer positions
		mset_surf.identify_monolayers(director)
		surf_select = mset_surf.universe.selectAtoms(selector)
		monoz = []
		for fr in whichframes:
			mset_surf.gotoframe(fr)
			surfpts = surf_select.coordinates()
			monoz.append([mean(surfpts[mset_surf.monolayer_residues[m],2],axis=0) for m in range(2)])
		#---midz and monoz provide the average midplane and average monolayer positions for each frame
		monoz = array(monoz)
		midz = mean(monoz,axis=1)
		
#---example for doing the coordinate shift and the binning in the z-dimension
if 'compute_z' in routine:
	aname = analysis_names[0]
	#---shift coordinate systems so the frame-wise midplane is at zero
	relpos = bilayer_shift()
	#---return a list of bin edges and the bin index for the bin flanking the monolayer
	#---choose 'fixed' for a fixed bin width or 'outside-inside' to get flush bins that might not be even
	if 0:
		binedges,monobins = binlister('fixed',mset=mset_surf,binw=5)
		binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-6,-5,10,20,30,32,40])
		binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-5,10])
		binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-5,0,15])
		binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-10,16])
		binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=list(arange(-5,40+1,5)))	
		binedges,monobins = binlister('fixed',mset=mset_surf,binw=5)
	binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-10,16])
	disctraj = discretizer_z(binedges)
	if 0:
		#---check the distribution
		plot_ion_distribution(disctraj,binedges)
		#---plot transition matrix
		plot_ion_transition_matrix(disctraj,binedges)
		#---plot residence times
		plot_ion_residence(disctraj,binedges,ignorezero=True,mset=mset_surf)
		plot_ion_residence(disctraj,binedges,ignorezero=True,mset=mset_surf,cumulative=False)
		plot_ion_time_correlations(disctraj,binedges,mset=mset_surf)

#---example for doing the coordinate shift and the binning relative to key phospholipids
if 'compute_radial_binary' in routine:
	aname = analysis_names[0]
	#---define selections
	#---Nb the following definitions need to be moved to the header, or at least get ion name from dictionary
	ionname = 'NA'
	buffer_size = 4
	selstring_refs = 'name O4 or name OP42 or name OP43 or name OP44 or name P4 or name O5 '+\
		'or name OP52 or name OP53 or name P5 or name O11 or name O12 or name O13 or name O14 or name P'
	selstring = 'name '+ionname+' and around '+str(buffer_size)+' ('+selstring_refs+')'
	expt_name = 'buffer'+str(buffer_size)+'-allheadgroup'
	#---discretize the trajectory
	frameslice = slice(0,400) #---set to None if you want the whole trajectory
	frameslice = None
	disctraj = discretizer_radial_binary(selstring,ionname,mset=mset,frameslice=frameslice)
	if 0:
		#---plot the residence time distributions
		plot_ion_residence_radial_binary(disctraj,'radial_binary',scale_by_time=False)
		#---new method that averages the distribution of leaving times for each t0
		plot_ion_residence_radial_binary_dev(disctraj,'radial_binary',normed=False)
		#---PDF of durations in the zones
		plot_ion_residence_radial_binary(disctraj,'radial_binary',scale_by_time=False)
		plot_ion_residence_radial_binary(disctraj,'radial_binary',scale_by_time=False,cumulative=True)
		plot_ion_residence_radial_binary_dev(disctraj,'radial_binary',normed=False,residence_exact=True)

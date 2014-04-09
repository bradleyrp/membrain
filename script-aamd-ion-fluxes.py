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
		'sysname_lookup':'membrane-v509-ions',
		'trajsel':'s6-kraken-md.part0018.40000-90000-10.ions.xtc',
		'sysname_lookup_struct':'membrane-v509-atomP',
		'trajsel_struct':'s6-kraken-md.part0018.40000-90000-10.atomP.xtc',
		'ionname':'NA',
		'director':director_aamd_symmetric}}
analysis_names = [
	'v532-20000-58000-100',
	'v530-30000-100000-100',
	'v531-20000-62000-100',
	'v509-40000-90000-10'
	][-2:-1]
routine = ['load','compute',][1:]

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

#---following code is on the chopping block ...
if 'compute_old' in routine:
      for aname in analysis_names:
		desired_binsize = 1
		upper_binws = [int(round((mset_surf.vec(i)[2]-monoz[i][0])/desired_binsize)) for i in whichframes]
		mid_binws = [int(round((monoz[i][0]-monoz[i][1])/desired_binsize)) for i in whichframes]
		lower_binws = [int(round((monoz[i][1])/desired_binsize)) for i in whichframes]
		#---old code frame-wise, mistaken because number of bins not consistent
		if 0:
			binedges = array([np.concatenate((
				linspace(0,monoz[i][0],lower_binws[i])[:-1],
				linspace(monoz[i][0],monoz[i][1],mid_binws[i])[:-1],
				linspace(monoz[i][1],mset_surf.vec(i)[2],upper_binws[i]))) 
				for i in whichframes])
			badframes = list(where([ionspos[j][:,2].max()>mset_surf.vec(j)[2] for j in range(5000)])[0])	
			whichframes2 = [i for i in whichframes if i not in badframes]
			disctraj  = array([[
				where(ionspos[j][:,2][i]<binedges[j][1:])[0][0] 
				for i in range(len(ionspos[0]))] for j in whichframes2]).T
		#---pick a standard number of bins for consistent zones
		bin_nums = [int(round(mean(i))) for i in [lower_binws,mid_binws,upper_binws]]
		badframes = list(where([ionspos[j][:,2].max()>mset_surf.vec(j)[2] for j in range(5000)])[0])	
		whichframes2 = [i for i in whichframes if i not in badframes]
		binedges = array([np.concatenate((
			linspace(0,monoz[i][1],bin_nums[0])[:-1],
			linspace(monoz[i][1],monoz[i][0],bin_nums[1])[:-1],
			linspace(monoz[i][0],mset_surf.vec(i)[2],bin_nums[2]))) 
			for i in whichframes])
		disctraj  = array([[
			where(ionspos[j][:,2][i]<binedges[j][1:])[0][0] 
			for i in range(len(ionspos[0]))] for j in whichframes2]).T
		#---quick calculation of a transition matrix
		if 0:
			disctrajt = disctraj.T
			tmat = zeros((len(binedges[0])-1,len(binedges[0])-1))
			for fr in range(len(disctrajt)-1):
				for i in range(len(disctrajt[fr])):
					tmat[disctrajt[fr][i],disctrajt[fr+1][i]] += 1
			#---plot the transition matrix
			plt.imshow(tmat,interpolation='nearest',origin='lower');plt.show()
		#---bleeding edge code makes an ugly plot of residence distribution times
		bincount = len(binedges[0])-1
		if 0:
			nbins = 100
			jumptimes = [where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0] for i in range(len(disctraj))]
			jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
			maxresid = max([max(i) for i in jumpdurations])
			jumphists = [histogram(jumpdurations[i],range=(0,maxresid),bins=nbins)[0] for i in range(len(jumpdurations))]
			jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
			residences = [[] for i in range(bincount)]
			for k in range(len(jumpbins)):
				for i in range(len(jumpbins[k])):
					residences[jumpbins[k][i]].append(jumpdurations[k][i])
			maxresid = max([max(i) for i in residences if i != []])
			resdists = [histogram(residences[i],range=(0,maxresid),bins=nbins)[0] for i in range(len(residences)) if i != []]	
			#---accounting
			binlabels = [int(round(i)) for i in (mean(binedges,axis=0)[1:]+mean(binedges,axis=0)[:-1])/2.]
			maxz = binedges[0][-1]
			meanz = mean(mean(monoz,axis=0))
			#tmp = [(i-meanz)+maxz*((i-meanz)<0)-maxz*((i-meanz)>maxz) for i in binlabels]
			howfar = [i-meanz-maxz*((i-meanz)>maxz/2.) for i in binlabels]
			#---plot
			ax = plt.subplot(111)
			for r in range(len(resdists)):
				dat = resdists[r]
				if howfar[r] > 0.: c = 'b'
				else: c = 'r'
				ax.plot(dat,c=c,lw=2,label=binlabels[r],
					alpha=(1-abs(howfar[r])/max(abs(array(howfar)))/1.1))
			ax.set_xlim((1,100))
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_title('residence time distributions v509')
			#ax.legend(loc=1)
			plt.show()
		if 0:
			#---plot some single-frame distributions
			for i in range(0,shape(disctraj)[1],1):
			#for i in [1]:
				hist,edges = histogram(disctraj[:,i],range=(0,bincount),bins=bincount);
				mids = array(edges[1:]+edges[:-1])/2.
				plt.plot(mids,hist);
			plt.show()
		if 0:
			#---average the histograms
			hists = []
			meanbinedges = mean(binedges,axis=0)	
			monobins = [array([abs(mean(monoz,axis=0)[i]-meanbinedges[j])
				for j in range(len(meanbinedges))]).argmin() for i in range(2)][::-1]
			for i in range(0,shape(disctraj)[1],1):
				hist,edges = histogram(disctraj[:,i],range=(0,bincount),bins=bincount);
				hists.append(hist)
			mids = array(edges[1:]+edges[:-1])/2.
			bw = mean(meanbinedges[1:]-meanbinedges[:-1]) # bin width in Angstroms.
			# plt.axvline(x=bw*edges[monobins[0]],ymin=0.,ymax=1.)
			# plt.axvline(x=bw*edges[monobins[1]],ymin=0.,ymax=1.)
			plt.plot(bw*mids,mean(hists,axis=0),'o-');
			plt.show()
			# Where is the histogram zero? Where is the center of reflection?
			zeroes = where(mean(hists,axis=0)==0)
			middle_bin = int(mean(zeroes))
		
#---UNDER CONSTRUCTION
#---re-work the discrete binning function, fix badframes, make flexible zones, etc

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
	'''This function generates a list of bin limits for a particular system. Feed this to discretizer.'''
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

def discretizer(binedges,ionspos=None):
	'''Discretizes an ion trajectory given a set of positions and a list of bins.'''
	if ionspos == None: ionspos = globals()['ionspos']
	disctraj  = array([[where(relpos[j,i]<binedges[j][1:])[0][0] for i in range(len(relpos[j]))] 
		for j in range(len(binedges))]).T
	return disctraj
	
def plot_ion_distribution(disctraj,binedges,fig=None):
	'''Generic plotting function for ion distributions by height from midplane according to binedges.'''
	if fig == None: fig = plt.figure()
	gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.0)
	ax = fig.add_subplot(gs[0])
	counts,edges = histogram(disctraj.flatten(),range=(0,len(binedges[0])-1),bins=len(binedges[0])-1,normed=True)
	meanedges = mean(binedges,axis=0)
	mids = meanedges[:-1]
	ax.bar(mids,counts,width=meanedges[1:]-meanedges[:-1])
	ax.grid(True)
	plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
	plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
	ax.set_xlabel(r'$z(x,y)\,\mathrm{(nm)}$',fontsize=fsaxlabel)
	ax.set_ylabel('frequency',fontsize=fsaxlabel)
	plt.show()
	return counts,edges
	
def plot_ion_transition_matrix(disctraj,binedges,fig=None):
	'''Plot the transition matrix for the discrete ion trajectories.'''
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
	
def plot_ion_residence(disctraj,binedges,fig=None,ignorezero=True,maxtime=100,
	mset=None,dofit=False,cumulative=True):
	'''Plot residence time distributions.'''
	if fig == None: fig = plt.figure(figsize=(8,10))
	if mset == None: mset = globals()['mset']
	#---fineness of the time bins is currently hard-coded
	ntimes = len(disctraj[0])
	jumptimes = [concatenate(([0],where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0],[len(disctraj[0])])) 
		for i in range(len(disctraj))]
	jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
	#---needs fixed
	#jumpdurations = [i[i!=0] for i in jumpdurations]
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

#---example for doing the coordinate shift and the binning
if 'compute' in routine:
	aname = analysis_names[0]
	#---shift coordinate systems so the frame-wise midplane is at zero
	relpos = bilayer_shift()
	#---return a list of bin edges and the bin index for the bin flanking the monolayer
	#---choose 'fixed' for a fixed bin width or 'outside-inside' to get flush bins that might not be even
	#binedges,monobins = binlister('fixed',mset=mset_surf,binw=5)
	#binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-6,-5,10,20,30,32,40])
	#binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-5,10])
	#binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-5,0,15])
	binedges,monobins = binlister('custom',mset=mset_surf,binw=5,custom_posts=[-10,16])
	disctraj = discretizer(binedges)
	#---check the distribution
	plot_ion_distribution(disctraj,binedges)
	#---plot transition matrix
	plot_ion_transition_matrix(disctraj,binedges)
	#---plot residence times
	plot_ion_residence(disctraj,binedges,ignorezero=True,mset=mset_surf)
	

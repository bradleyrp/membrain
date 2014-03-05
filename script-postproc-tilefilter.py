#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

from scipy.spatial.distance import *

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'protein_traj':None,
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$'},
	'v700-500000-600000-200':
		{'sysname':'membrane-v700','sysname_lookup':None,
		'trajsel':'u1-lonestar-longrun/md.part0009.500000-700000-200.xtc',
		'protein_traj':None},
	'v701-60000-160000-200':
		{'sysname':'membrane-v701','sysname_lookup':None,
		'trajsel':'s8-lonestar/md.part0003.60000-160000-200.xtc',
		'protein_traj':None},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'protein_traj':None,
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$'},
	'v550-4000000-500000-160':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc',
		'protein_traj':None},
	'v550-300000-400000-200':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'protein_traj':'v614-120000-220000-200',
		'custom_protein_shifts':['unshifted','peak','valley'],
		'label':r'$\mathrm{control}$'}}

routine = ['topogareas','topogareas_singleplots'][-1:]
analysis_names = ['v614-120000-220000-200','v612-75000-175000-200','v550-300000-400000-200'][:]
plot_reord = analysis_names
bigname = 'v614-v612-v550'

#---parameters
binwid = 5.

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = []
	topogareas = []
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		msetfile = 'pkl.structures.'+specname_guess(sysname,trajsel)+'.pkl'
		msets.append(unpickle(pickles+msetfile))
		if protein_traj == None:
			topogarea = unpickle(pickles+'pkl.topography_transform.'+\
				specname_guess(sysname,trajsel)+'.pkl').data
		else:
			topogarea = []
			for shiftname in custom_protein_shifts:
				readdat = unpickle(pickles+'pkl.topography_transform.'+\
					specname_guess(sysname,trajsel)+'.'+shiftname+'.pkl')
				if readdat == None:	topogarea.append(None)
				else: topogarea.append(readdat.data)
		if topogarea == None or all([i==None for i in topogarea]):
			print 'status: generating topography_transform data'
			if protein_traj == None:
				mset = msets[-1]
				tilefiltdat = []
				for fr in range(len(mset.surf)):
					print fr
					protpts = mset.protein[fr][:,:2]
					surfpts = mset.wrappbc(mset.unzipgrid(array(mset.surf[fr]),vecs=mset.vec(0)),
						vecs=mset.vec(fr),mode='nine')
					cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2])
					rawdat = array([np.min(cd,axis=0),surfpts[:,2]]).T
					tilefiltdat.append(rawdat)
				print 'status: pickling the result'
				result_data = MembraneData('topography_transform')
				result_data.data = array(array(tilefiltdat))
				for i in analysis_descriptors[aname]: 
					result_data.addnote([i,(analysis_descriptors[aname])[i]])
				pickledump(result_data,'pkl.topography_transform.'+specname_guess(sysname,trajsel)+'.pkl',
					directory=pickles)
				topogareas.append(result_data.data)
			#---if control use different protein points
			else:
				topogareas_control = []
				for shiftname in custom_protein_shifts:
					msetprotfile = 'pkl.structures.'+specname_guess(\
						(analysis_descriptors[protein_traj])['sysname'],
						(analysis_descriptors[protein_traj])['trajsel'])+'.pkl'
					mset_prot = unpickle(pickles+msetprotfile)
					mset = msets[-1]
					tilefiltdat = []
					if shiftname == 'unshifted':
						shift = array([0.,0.])
					elif shiftname == 'peak':
						maxpos = unravel_index(mean(mset.surf,axis=0).argmax(),mset.surf[0].shape)
						maxxy = [maxpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
						protein_com = mean(mean(mset_prot.protein,axis=0),axis=0)[:2]
						shift = maxxy - protein_com
						print shift
					elif shiftname == 'valley':
						minpos = unravel_index(mean(mset.surf,axis=0).argmin(),mset.surf[0].shape)
						minxy = [minpos[i]*mean(mset.vecs,axis=0)[i]/mset.griddims[i] for i in range(2)]
						protein_com = mean(mean(mset_prot.protein,axis=0),axis=0)[:2]
						shift = minxy - protein_com
						print shift
					for fr in range(len(mset.surf)):
						print fr
						protpts = mset_prot.protein[fr][:,:2]
						surfpts = mset.wrappbc(mset.unzipgrid(array(mset.surf[fr]),vecs=mset.vec(0)),
							vecs=mset.vec(fr),mode='nine')
						cd = scipy.spatial.distance.cdist(protpts,surfpts[:,:2]+shift)
						rawdat = array([np.min(cd,axis=0),surfpts[:,2]]).T
						tilefiltdat.append(rawdat)
					print 'status: pickling the result'
					result_data = MembraneData('topography_transform')
					result_data.data = array(array(tilefiltdat))
					for i in analysis_descriptors[aname]: 
						result_data.addnote([i,(analysis_descriptors[aname])[i]])
					pickledump(result_data,'pkl.topography_transform.'+\
						specname_guess(sysname,trajsel)+'.'+shiftname+'.pkl',
						directory=pickles)
					topogareas_control.append(result_data.data)
				topogareas.append(topogareas_control)
		else:
			topogareas.append(topogarea)
			
#---compute tilefiltered areas
if 'topogareas_singleplots' in routine:
	pnum = 0
	npanels = sum([(1 if i == None else 3)
		for i in [(analysis_descriptors[aname])['protein_traj'] for aname in analysis_names]])
	if 'all_areacurves' not in globals(): all_areacurves = [[] for i in range(npanels)]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,npanels)
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		m = analysis_names.index(aname)
		print aname
		mset = msets[m]
		if protein_traj == None: datlist = [topogareas[m]]
		else: datlist = topogareas[m]
		for p in range(len(datlist)):
			topogarea = datlist[p]
			cutoff = min(mean(mset.vecs,axis=0)[:2])
			#---training wheels to accelerate the plot-tweaking
			if all_areacurves[pnum] == []:
				areacurves = []
				for rawdat in topogarea:
					areacurves.append([float(sum(rawdat[rawdat[:,0]<i,1]>0.))/sum(rawdat[:,0]<i) 
						for i in arange(binwid,cutoff,binwid)])
				areacurves = array(areacurves)
				all_areacurves[pnum] = areacurves
			else:
				areacurves = all_areacurves[pnum]
			ax = fig.add_subplot(gs[pnum])
			ax.axhline(y=0.5,linewidth=2,color='k')
			ax.fill_between(arange(binwid,cutoff,binwid)/10.,
				mean(areacurves,axis=0)-std(areacurves,axis=0),
				mean(areacurves,axis=0)+std(areacurves,axis=0),facecolor='k',alpha=0.5)
			ax.plot(arange(binwid,cutoff,binwid)/10.,mean(areacurves,axis=0),'k-',lw=3,
				markersize=2)
			ax.grid(True)
			ax.set_xlim((0,cutoff/10.))
			ax.set_ylim((0,1.0))
			#---plot details
			ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=3,prune='both'))
			ax.set_title(label,fontsize=fsaxlabel)
			plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
			plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
			ax.set_xlabel(r'$r_{min}\:(\mathrm{nm})$',fontsize=fsaxlabel)
			if pnum == 0: ax.set_ylabel(r'$A_{z>0}/A$',fontsize=fsaxlabel)
			if pnum > 0: ax.set_yticklabels([])
			pnum += 1
	plt.savefig(pickles+'fig-tilefilter-areas-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.clf()	
	plt.cla()
	plt.close()


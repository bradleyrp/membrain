#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---plan
analysis_descriptors = {
	'v614-120000-220000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'trajsel':'s9-lonestar/md.part0004.120000-220000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4\,(v2)}$',
		'nprots':4,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'topogcorr_pkl':'pkl.topogcorr.membrane-v614-s9-lonestar-120000-220000-200.pkl'},
	'v612-75000-175000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'t4-lonestar/md.part0007.75000-175000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1\,(v2)}$',
		'nprots':1,
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'whichframes':slice(None,None)},
	'v614-40000-140000-200':
		{'sysname':'membrane-v614','sysname_lookup':None,
		'trajsel':'s6-sim-lonestar/md.part0002.40000-140000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}4}$',
		'nprots':4,
		'whichframes':slice(None,None),
		'protein_pkl':None,
		'custom_topogcorr_specs':None},
	'v612-10000-80000-200':
		{'sysname':'membrane-v612','sysname_lookup':None,
		'trajsel':'s9-trestles/md.part0003.10000-80000-200.xtc',
		'label':r'$\mathrm{{ENTH}\ensuremath{\times}1}$',
		'nprots':1,
		'protein_pkl':None,
		'custom_topogcorr_specs':None,
		'whichframes':slice(None,None)},		
	'v550-300000-400000-200':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'s0-trajectory-full/md.part0006.300000-400000-200.xtc',
		'label':r'$\mathrm{control}$',
		'nprots':1,
		'whichframes':slice(0,500),
		'protein_pkl':'pkl.structures.membrane-v612.t4-lonestar.md.part0007.75000-175000-200.pkl',
		'custom_topogcorr_specs':None,
		'custom_protein_shifts':['peak','valley']},
	'v550-400000-500000-160':
		{'sysname':'membrane-v550','sysname_lookup':None,
		'trajsel':'v1-lonestar/md.part0010.400000-500000-160.xtc',
		'label':r'$\mathrm{control\,(v2)}$',
		'nprots':1,
		'whichframes':slice(0,500),
		'protein_pkl':'pkl.structures.membrane-v612.t4-lonestar.md.part0007.75000-175000-200.pkl',
		'custom_topogcorr_specs':None,
		'custom_protein_shifts':['peak','valley']}}
		
#---analysis menu
analysis_names = ['v614-120000-220000-200','v612-75000-175000-200','v550-400000-500000-160',
	'v614-40000-140000-200','v612-10000-80000-200','v550-300000-400000-200'][:]
plot_reord = analysis_names
routine = ['plot','plot_mean_z','plot_gaussian'][:]
bigname = 'v614-v612-v550-old-new'
panel_nrows = 2

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def curvcalc(z,lenscale):
	'''Calculate mean and Gaussian curvature directly.'''
	zy, zx  = numpy.gradient(z,lenscale)
	zxy, zxx = numpy.gradient(zx,lenscale)
	zyy, _ = numpy.gradient(zy,lenscale)
	H = (zx**2 + 1)*zyy - 2*zx*zy*zxy + (zy**2 + 1)*zxx
	H = -H/(2*(zx**2 + zy**2 + 1)**(1.5))
	K = ((zxx*zyy)-(zxy)**2)
	K = -K/(2*(zx**2 + zy**2 + 1)**(1.5))
	return [H,K]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		msets.append(unpickle(pickles+'pkl.structures.'+specname_guess(sysname,trajsel)+'.pkl'))

#---compute curvatures
if 'calc' in routine or 'curvsm_all' not in globals():
	curvsm_all = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		curvsm = []
		for i in range(len(mset.surf)):
			lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
			curvsm.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[0])
		curvsm_all.append(array(curvsm))
	curvsk_all = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		curvsk = []
		for i in range(len(mset.surf)):
			lenscale = mean(mset.vecs,axis=0)[0]/10./(shape(mset.surf[0])[0]+1)
			curvsk.append(curvcalc(list(array(mset.surf[i]).T),lenscale)[1])
		curvsk_all.append(array(curvsk))


#---plot curvatures
if 'plot' in routine:
	extremum = max([abs(mean(i,axis=0).max()) for i in curvsm_all]+[abs(mean(i,axis=0).min()) 
		for i in curvsm_all])
	fig = plt.figure()
	if panel_nrows > 1:
		gs = gridspec.GridSpec(panel_nrows,int(ceil(float(len(msets))/panel_nrows)))
	else:
		gs = gridspec.GridSpec(1,len(msets))
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None: mset_protein = unpickle(pickles+protein_pkl)
		else: mset_protein = mset
		if panel_nrows > 1:
			ax = fig.add_subplot(gs[m/(len(msets)/panel_nrows),int(m%(floor(len(msets)/panel_nrows)))])
		else:
			ax = fig.add_subplot(gs[m])
		ax.set_title(label,fontsize=fsaxtitle)
		im = plotter2d(ax,mset,dat=mean(curvsm_all[m],axis=0),tickshow=True,
			label_style='xy',lims=[-1*extremum,extremum],lognorm=False,
			cmap=mpl.cm.RdBu_r,fs=fsaxlabel)
		if protein_pkl == None:
			print a
			plothull(ax,mean(mset_protein.protein,axis=0),mset=mset_protein,c='k',subdivide=nprots)
		if m != 0: ax.set_ylabel('')
		if m == len(msets)-1 or (panel_nrows > 1 and 
			int(m%(floor(len(msets)/panel_nrows))) == (len(msets)/panel_nrows-1)):
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins.set_ylabel(r'$\left\langle H(x,y)\right\rangle \:(\mathrm{{nm}^{-1}})$',
				fontsize=fsaxlabel,rotation=270)
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.2,w_pad=0.2)
	plt.savefig(pickles+'fig-curve_map-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.show()

#---plot average structure virtually identical to the codeblock above
if 'plot_mean_z' in routine:
	extremum = max([abs(mean(i,axis=0).max()) for i in [msets[j].surf for j in range(len(msets))]]+\
		[abs(mean(i,axis=0).min()) for i in [msets[j].surf for j in range(len(msets))]])
	fig = plt.figure()
	if panel_nrows > 1:
		gs = gridspec.GridSpec(panel_nrows,int(ceil(float(len(msets))/panel_nrows)))
	else:
		gs = gridspec.GridSpec(1,len(msets))
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None: mset_protein = unpickle(pickles+protein_pkl)
		else: mset_protein = mset
		if panel_nrows > 1:
			ax = fig.add_subplot(gs[m/(len(msets)/panel_nrows),int(m%(floor(len(msets)/panel_nrows)))])
		else:
			ax = fig.add_subplot(gs[m])
		ax.set_title(label,fontsize=fsaxtitle)
		im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0),tickshow=True,
			label_style='xy',lims=[-1*extremum,extremum],lognorm=False,
			cmap=mpl.cm.RdBu_r,fs=fsaxlabel)
		if protein_pkl == None:
			print a
			plothull(ax,mean(mset_protein.protein,axis=0),mset=mset_protein,c='k',subdivide=nprots)
		if m != 0: ax.set_ylabel('')
		if m == len(msets)-1 or (panel_nrows > 1 and 
			int(m%(floor(len(msets)/panel_nrows))) == (len(msets)/panel_nrows-1)):
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins.set_ylabel(r'$\left\langle H(x,y)\right\rangle \:(\mathrm{{nm}^{-1}})$',
				fontsize=fsaxlabel,rotation=270)
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.2,w_pad=0.2)
	plt.savefig(pickles+'fig-mean_z_map-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.show()

#---plot curvatures
if 'plot_gaussian' in routine:
	extremum = max([abs(mean(i,axis=0)[1:-1,1:-1].max()) for i in curvsk_all]+\
		[abs(mean(i,axis=0)[1:-1,1:-1].min()) for i in curvsk_all])
	fig = plt.figure()
	if panel_nrows > 1:
		gs = gridspec.GridSpec(panel_nrows,int(ceil(float(len(msets))/panel_nrows)))
	else:
		gs = gridspec.GridSpec(1,len(msets))
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		m = analysis_names.index(a)
		mset = msets[m]
		if protein_pkl != None: mset_protein = unpickle(pickles+protein_pkl)
		else: mset_protein = mset
		if panel_nrows > 1:
			ax = fig.add_subplot(gs[m/(len(msets)/panel_nrows),int(m%(floor(len(msets)/panel_nrows)))])
		else:
			ax = fig.add_subplot(gs[m])
		ax.set_title(label,fontsize=fsaxtitle)
		im = plotter2d(ax,mset,dat=mean(curvsk_all[m],axis=0)[1:-1,1:-1],tickshow=True,
			label_style='xy',lims=[-1*extremum,extremum],lognorm=False,
			cmap=mpl.cm.RdBu_r,fs=fsaxlabel)
		if protein_pkl == None:
			print a
			plothull(ax,mean(mset_protein.protein,axis=0),mset=mset_protein,c='k',subdivide=nprots)
		if m != 0: ax.set_ylabel('')
		if m == len(msets)-1 or (panel_nrows > 1 and 
			int(m%(floor(len(msets)/panel_nrows))) == (len(msets)/panel_nrows-1)):
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
			axins.set_ylabel(r'$\left\langle K(x,y)\right\rangle \:(\mathrm{{nm}^{-2}})$',
				fontsize=fsaxlabel,rotation=270)
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.2,w_pad=0.2)
	plt.savefig(pickles+'fig-curve_map_gaussian-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.show()












		

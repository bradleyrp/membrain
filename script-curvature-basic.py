#!/usr/bin/python -i

interact = True
from membrainrunner import *
execfile('locations.py')

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---get dictionaries
execfile('header-cgmd.py')
	
#---analysis menu
analysis_names = [
	'v614-120000-220000-200',
	'v612-75000-175000-200',
	'v550-400000-500000-160',
	'v614-40000-140000-200',
	'v612-10000-80000-200',
	'v550-300000-400000-200',
	'v616-210000-310000-200',
	'v701-60000-160000-200',
	'v550-300000-400000-200',
	][-3:]
plot_reord = analysis_names
routine = [
	'plot',
	'plot_mean_z',
	'plot_gaussian',
	'video_height',
	'video_height_cat',
	'extrema_dynamics',
	][1:2]
bigname = 'v616-v614-v612-v550'
bigname = 'v616-combine'
bigname = 'v616-v701-v550'
panel_nrows = [1,2][0]

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
	
def gridprop_panel_plot(dat,fig,fr=None,cmap=None,vmax=None,vmin=None,
	altdat=None,extrarow=False,nine=False,args=None,fill=None,sidebar=True,protalpha=1,nrows=None):
	'''Function which plots a map with proteins.'''
	#---process pass-through argument list from plotmov if available
	if args != None:
		for arg in args:
			if arg[0] == 'nine': nine = arg[1]
	#---settings
	panels = len(dat)
	if cmap == None: cmap = mpl.cm.RdBu_r
	if fr == None: fr = 0
	if nrows == None: nrows = 1
	#---axes
	axeslist = []
	#gs = gridspec.GridSpec((2 if extrarow else 1),panels,wspace=0.0,hspace=0.0)
	#---prepare rows and columns
	ncols = int(ceil(float(panels)/nrows))
	gs = gridspec.GridSpec(nrows,ncols,wspace=0.0,hspace=0.0)
	poslist = [[i,j] for i in range(nrows) for j in range(ncols)]
	for p in range(panels):
		#if extrarow: ax = fig.add_subplot(gs[1,p])
		#else: ax = fig.add_subplot(gs[p])
		ax = fig.add_subplot(gs[poslist[p][0],poslist[p][1]])
		axeslist.append(ax)
		im = plotter2d(ax,msets[p],dat=dat[p],tickshow=True,
			label_style='xy',lims=[vmin,vmax],lognorm=False,
			cmap=cmap,fs=fsaxlabel,nine=nine,tickskip=(25 if nine else 10))
		if altdat != None and altdat[p].protein != []:
			if nine:
				vecs=mean(msets[p].vecs,axis=0)
				for shift in [[i*vecs[0],j*vecs[1],0] for i in range(3) for j in range(3)]:
					plothull(ax,msets[p].protein[fr]+shift,griddims=msets[0].griddims,
						vecs=mean(msets[p].vecs,axis=0),subdivide=nnprots[p],alpha=protalpha,
						c='k',fill=fill)
			else:	
				plothull(ax,msets[p].protein[fr],griddims=msets[0].griddims,
					vecs=mean(msets[p].vecs,axis=0),subdivide=nnprots[p],alpha=protalpha,
					c='k',fill=fill)
		#---if plotting the nine cells, draw a box in the middle
		if nine:
			m,n = msets[p].griddims
			boxlims = [[m,2*m],[n,2*n]]
			ax.axvline(x=boxlims[0][0]-1,ymin=1./3,ymax=2./3,linewidth=1,c='k')
			ax.axvline(x=boxlims[0][1],ymin=1./3,ymax=2./3,linewidth=1,c='k')
			ax.axhline(y=boxlims[1][0]-1,xmin=1./3,xmax=2./3,linewidth=1,c='k')
			ax.axhline(y=boxlims[1][1],xmin=1./3,xmax=2./3,linewidth=1,c='k')
		#---more labels
		ax.set_xlabel(r'$x\:(\mathrm{nm})$',fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		#---added the following for flush plots
		if nrows == 1:
			if p > 0: ax.set_yticklabels([])
			if p == 0: ax.set_ylabel(r'$y\:(\mathrm{nm})$',fontsize=fsaxlabel)
			else: ax.set_ylabel('')
		else:
			if poslist[p][1] > 0: 
				ax.set_ylabel('')
				ax.set_yticklabels([])
			if poslist[p][0] < ncols-1: 
				ax.set_xlabel('')
				ax.set_xticklabels([])
		if 0: ax.set_title(labels[p],fontsize=fsaxtitle)
	if sidebar:
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.,0.,1.,1.),
			bbox_transform=ax.transAxes,
			borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(axins.get_xticklabels(),fontsize=fsaxlabel)
		axins.set_ylabel(r'$\mathsf{z(x,y)\:(nm)}$',fontsize=fsaxlabel,rotation=270)
	return [gs,axeslist]
	
#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load
if 'msets' not in globals():
	msets = []
	nnprots = []
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		msets.append(unpickle(pickles+'pkl.structures.'+spacetag+specname_guess(sysname,trajsel)+'.pkl'))
		nnprots.append(nprots)

#---make cool videos of the height
if 'video_height' in routine:
	#---standard video routine modified from the stressmap code
	#---find extrema to normalize the color scale between panels
	extremum = max(max([array(i.surf).max() for i in msets]),
		max([abs(array(i.surf).min()) for i in msets]))/10.
	#---test the plotting function
	gridprop_panel_plot([array(msets[i].surf[0])/10. for i in range(len(msets))],
		plt.figure(),vmin=-extremum,vmax=extremum,fr=0,altdat=msets,nine=True,nrows=panel_nrows)
	plt.show()
	#---prepare the data set
	vidpanels = [[array(msets[i].surf[j])/10. 
		for j in range(len(msets[0].surf))] 
		for i in range(len(msets))]
	#---generate the video
	#---previously used figsize 12,6 with nrows = 1 for a 1x4 panel plot
	#---previously used figsize 12,6 with nrows = 2 for a 2x2 panel plot, with a slight gap between rows
	plotmov(vidpanels,'heightmap-'+bigname,altdat=msets,panels=len(msets),
		plotfunc='gridprop_panel_plot',whitezero=True,lims=[-extremum,extremum],
		args=[['nine',True]],slowdown=2,figsize=(10,10),lowres=True,nrows=panel_nrows)
		
#---make cool videos of the height
if 'video_height_cat' in routine:
	#---standard video routine modified from the stressmap code
	#---find extrema to normalize the color scale between panels
	extremum = max(max([array(i.surf).max() for i in msets]),
		max([abs(array(i.surf).min()) for i in msets]))/10.
	#---test the plotting function
	msetcat = msets[0]
	msetcat.surf = msetcat.surf + msets[1].surf
	msetcat.protein = msetcat.protein + msets[1].protein
	msetcat.vecs = msetcat.vecs + msets[1].vecs
	msetcat.vecs_index = list(array(msets[1].vecs_index)+len(msetcat.vecs))
	gridprop_panel_plot([array(msetcat.surf[0])/10.],
		plt.figure(),vmin=-extremum,vmax=extremum,fr=0,altdat=msets,nine=True,nrows=panel_nrows)
	plt.show()
	#---prepare the data set
	vidpanels = [[array(msetcat.surf[j])/10. 
		for j in range(len(msets[0].surf))]]
	#---generate the video
	#---previously used figsize 12,6 with nrows = 1 for a 1x4 panel plot
	#---previously used figsize 12,6 with nrows = 2 for a 2x2 panel plot, with a slight gap between rows
	plotmov(vidpanels,'heightmap-'+bigname,altdat=[msetcat],panels=1,
		plotfunc='gridprop_panel_plot',whitezero=True,lims=[-extremum,extremum],
		args=[['nine',True]],slowdown=2,figsize=(10,10),lowres=True,nrows=panel_nrows)

#---compute curvatures
if 'extrema_dynamics' in routine:
	fig = plt.figure(figsize=(12,6))
	extrema_pts = []
	for mset in msets:
		m,n = mset.griddims
		minpts = zeros((m,n))
		maxpts = zeros((m,n))
		for fr in range(len(mset.surf)):
			xi,yi = unravel_index(mset.surf[fr].argmin(),mset.surf[fr].shape)
			minpts[xi][yi] += 1
			xi,yi = unravel_index(mset.surf[fr].argmax(),mset.surf[fr].shape)
			maxpts[xi][yi] += 1
		if 0: extrema_pts.append(maxpts/sum(maxpts))
		extrema_pts.append((maxpts/mean(sum(maxpts)+sum(minpts))-minpts/mean(sum(maxpts)+sum(minpts))))
	topval = max([max([i.max(),abs(i.min())]) for i in extrema_pts])/4.
	gs,axlist = gridprop_panel_plot(extrema_pts,plt.figure(),vmin=-topval,vmax=topval,cmap=mpl.cm.RdBu_r,
		fr=0,altdat=msets,nine=True,sidebar=False,fill=False,protalpha=0.5)
	plt.savefig(pickles+'fig-peak_valley_positions-'+bigname+'.png',dpi=300,bbox_inches='tight')
	plt.show()

#---compute curvatures
if 'calc' in routine or ('curvsm_all' not in globals() and 
	('plot_mean_curvature' in routine or 'plot_gaussian_curvature' in routine)):
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
if 'plot_mean_curvature' in routine:
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
		#if protein_pkl != None: mset_protein = unpickle(pickles+protein_pkl)
		#if nprots == 0:
		#else: mset_protein = mset
		mset_protein = mset
		if panel_nrows > 1:
			ax = fig.add_subplot(gs[m/(len(msets)/panel_nrows),int(m%(floor(len(msets)/panel_nrows)))])
		else:
			ax = fig.add_subplot(gs[m])
		ax.set_title(label,fontsize=fsaxtitle-2)
		im = plotter2d(ax,mset,dat=mean(mset.surf,axis=0),tickshow=True,
			label_style='xy',lims=[-1*extremum,extremum],lognorm=False,
			cmap=mpl.cm.RdBu_r,fs=fsaxlabel)
		#if protein_pkl == None:
		#	print a
		#	plothull(ax,mean(mset_protein.protein,axis=0),mset=mset_protein,c='k',subdivide=nprots)
		if nprots > 0:
			plothull(ax,mean(mset_protein.protein,axis=0),mset=mset_protein,c='k',subdivide=nprots)
		if m != 0: ax.set_ylabel('')
		if m == len(msets)-1 or (panel_nrows > 1 and 
			int(m%(floor(len(msets)/panel_nrows))) == (len(msets)/panel_nrows-1)):
			axins = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.,0.,1.,1.),
				bbox_transform=ax.transAxes,
				borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			plt.setp(axins.get_yticklabels(),fontsize=fsaxlabel-4)
			axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle \:(\mathrm{{nm}})$',
				fontsize=fsaxlabel-2,rotation=270)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel-6)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel-6)
	fig.set_size_inches(fig.get_size_inches()[0]*1.5,fig.get_size_inches()[1]*1.5)
	gs.tight_layout(fig,h_pad=0.2,w_pad=0.2)
	plt.savefig(pickles+'fig-mean_z_map-'+bigname+'.png',dpi=500,bbox_inches='tight')
	plt.show()

#---plot curvatures
if 'plot_gaussian_curvature' in routine:
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












		

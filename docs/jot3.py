#!/usr/bin/python
showplots = True
def plotter_undulation_residuals(thresholds,comp,comp_cgmd,toptitle,a0,filename_descriptor):
	'''Plot a summary of the undulation residuals across systems.'''
	fig = plt.figure(figsize=((4,3*len(thresholds)) if len(thresholds)>1 else (8,8)))
	gs = gridspec.GridSpec(len(thresholds),2,wspace=0.1,hspace=0.1,
		width_ratios=[len(meso_c0s),len(cgmd_list)])
	cmap = mpl.cm.jet
	#---plots loop over rows where each row has the comparison for a different wavevector cutoff
	for thresh in thresholds:
		residual_comparisons = comp[thresholds.index(thresh)]
		residual_comparisons_cgmd = comp_cgmd[thresholds.index(thresh)]
		axl = plt.subplot(gs[thresholds.index(thresh),0])
		axr = fig.add_subplot(gs[thresholds.index(thresh),1])
		max_resid = max([residual_comparisons_cgmd.max(),residual_comparisons_cgmd.max()])
		im = axl.imshow(residual_comparisons,interpolation='nearest',cmap=cmap,
			vmax=max_resid,vmin=0)
		axl.set_xticks(range(len(meso_c0s)))
		axl.set_yticks(range(len(meso_c0s)))
		axl.set_xticklabels(['{:.3f}'.format(i/a0) for i in meso_c0s])
		axl.set_yticklabels(['{:.3f}'.format(i/a0) for i in meso_c0s])
		for label in im.axes.xaxis.get_ticklabels():
			label.set_rotation(90)
		axl.set_xlim((0,len(meso_c0s)))
		im = axr.imshow(residual_comparisons_cgmd,interpolation='nearest',cmap=cmap,
			vmax=max_resid,vmin=0)
		tmp = array(comp_cgmd)[0].T
		for k in range(4):
			for i in range(len(tmp)):
				pos = argsort(tmp[i])[k]
				axr.text(i,pos,str(k),fontsize=16,horizontalalignment='center',va='center',weight='bold')
		
		axr.set_xticks(range(len(cgmd_list)))
		axr.set_xticklabels([(analysis_descriptors[name])['label'] for name in cgmd_list])
		for label in im.axes.xaxis.get_ticklabels():
			label.set_rotation(90)
		plt.setp(axr.get_yticklabels(), visible=False)
		axins = inset_axes(axr,width=1./len(cgmd_list)/2.,height="100%",loc=3,
			bbox_to_anchor=(1.2,0.0,1.,1.),
			bbox_transform=axr.transAxes,
			borderpad=0)
		boundaries = linspace(0,max_resid,21)
		ticks = [float('{:.3f}'.format(i)) for i in boundaries]
		norm = mpl.colors.BoundaryNorm(ticks,cmap.N)
		cbar = mpl.colorbar.ColorbarBase(axins,cmap=cmap,ticks=ticks,boundaries=boundaries,
			norm=norm,spacing='proportional',format='%03f')
		cbar.set_ticklabels([str(i) for i in ticks])
		axl.set_xlim((-0.5,len(meso_c0s)-1+0.5))
		axr.set_xlim((-0.5,len(cgmd_list)-1+0.5))
		axl.set_title(toptitle,fontsize=fsaxtitle)
		axl.set_ylabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
		axl.set_xlabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
	#---save
	if 0:
		plt.savefig(pickles+'fig-bilayer-couple-meta-'+filename_descriptor+\
			'-'.join([i.split('-')[0] for i in cgmd_avail]+meso_avail)+\
			'.png',
			bbox_inches='tight',dpi=300)
	if showplots: plt.show()
	plt.close(fig)
sn = 0
plotter_undulation_residuals(thresholds,comp[sn],comp_cgmd[sn],toptitles[sn],a0,fnames[sn]+('shifted-' if shift_curves else ''))

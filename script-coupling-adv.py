#!/usr/bin/python

#---allocate if empty
if 'msets' not in globals(): msets = []
if 'mscs' not in globals(): mscs = []
if 'collect_c0s' not in globals(): collect_c0s = []

#---load and interpolate

lenscale = 1.0
for a in analysis_names:
	for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
	if 'mset' in globals(): del mset
	mset = MembraneSet()
	if simtype == 'meso':
		params = analysis_descriptors[a]
		#---? type on rundir is int but returns as a float
		pklname = 'pkl.structures.meso.'+\
			params['callsign']+'-'+\
			'C_0-'+str(params['C_0'])+'-'+\
			'len-'+str(params['lenscale'])+'-'+\
			'rundir-'+str(int(params['rundir']))+\
			'.pkl'
		print 'unpickle '+str(pklname)
		mset = unpickle(pickles+pklname)
		print array(mset.surf[0]).max()
		c0s = mset.getdata('c0map').data
		collect_c0s.append(c0s)
		msets.append(mset)
	elif simtype == 'md':
		status('status: loading from MD')
		if 'mset' in globals(): del mset
		mset = unpickle(pickles+locate)
		msets.append(mset)
		#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
		c0ask = (analysis_descriptors[analysis_names[1]])['C_0']
		hypo[0] = c0ask
		#---compute hypothetical curvature field for the MD simulation
		if hypo != None:
			vecs = mean(mset.vecs,axis=0)
			m,n = mset.griddims
			#---getgrid is xyz points in nm
			getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
				for i in linspace(0,vecs[0]/mset.lenscale,m)])
			#---convert everything to nm
			#---recall z0,c0,x0,y0,sx,sy,th = params
			#---params sets C0 in native units, x0,y0 in proportional units, and sx,sy in nm
			params = [0,hypo[0],
				vecs[0]*hypo[1]/mset.lenscale,
				vecs[1]*hypo[2]/mset.lenscale,
				hypo[3],hypo[4],
				hypo[5]]
			c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)]
				for i in range(m)])
			collect_c0s.append([c0hypo for i in range(len(mset.surf))])
		else:
			collect_c0s.append([])
	elif simtype == 'meso_precomp':
		if 'mset' in globals(): del mset
		mset = unpickle(pickles+locate)
		#---for precomp simulations these units are relative to a0 and already in a reglar grid
		collect_c0s.append(mset.getdata('c0map').data)
		msets.append(mset)

#---calculate mode couplings

#---match mesoscale length scales to an MD simulation
if match_scales != None:
	ref_ind = analysis_names.index(match_scales[0])
	move_ind = analysis_names.index(match_scales[1])
	#---matching the average box vectors here by average and not maximum
	lenscale = mean(mean(msets[move_ind].vecs,axis=0)[:2])/\
		(mean(mean(msets[ref_ind].vecs,axis=0)[:2])/msets[ref_ind].lenscale)
	for a in analysis_names:
		if (analysis_descriptors[a])['simtype'] == 'meso' or \
			(analysis_descriptors[a])['simtype'] == 'meso_precomp':
			#analysis_descriptors[a]['lenscale'] = lenscale
			for i in analysis_descriptors[a]: 
				if i != 'lenscale': vars()[i] = (analysis_descriptors[a])[i]
			msets[analysis_names.index(a)].lenscale = lenscale
	#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
	hypo[0] = c0ask
	#---reset the hypothetical C0 field according to the new scaling (replaces the calculation above)
	for a in analysis_names:
		for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
		if analysis_descriptors[a]['simtype'] == 'md':
			anum = analysis_names.index(a)
			mset = msets[anum]
			vecs = mean(mset.vecs,axis=0)
			m,n = mset.griddims
			#---getgrid is xyz points in nm
			getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
				for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
			#---needs checked, "key step used to have a 0.5*hypo[0] here possibly due to convention"
			params = [0,
				hypo[0]*msets[1].lenscale/mset.lenscale,
				vecs[0]*(1+hypo[1])/mset.lenscale,
				vecs[1]*(1+hypo[2])/mset.lenscale,
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				sqrt(r_2)/msets[1].lenscale/sqrt(2),
				hypo[5]]
			#---handle curvature fields at the mesoscale which cross the PBC boundary
			c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
				for j in range(3*n)] for i in range(3*m)])
			c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
				for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
			collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]

#---calculate coupled modes
for a in analysis_names:
	for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
	m = analysis_names.index(a)
	if 'msc' in globals(): del msc
	msc = ModeCouple()
	msc.calculate_mode_coupling(msets[m],collect_c0s[m])
	mscs.append(msc)
hypo = (analysis_descriptors[batch_cgmd])['hypo']
hypo[0] = c0ask
if 'masterplot' not in routine and 'simple_summary' in routine: spectrum_summary()
if 'batch_override' in globals() and batch_override:
	if 0:
		if '-'.join(analysis_names) not in master_spectrum_dict.keys():
			master_spectrum_dict['-'.join(analysis_names)] = {\
				'c0ask':c0ask,
				'cgmd_qs':mscs[0].qmagst,
				'cgmd_t2d':mscs[0].t2d,
				'lenscale':msets[0].lenscale,
				}
	for s in range(2):
		if analysis_names[s] not in master_spectrum_dict.keys():
			master_spectrum_dict[analysis_names[s]] = {\
				'c0ask':c0ask,
				('cgmd_qs' if analysis_descriptors[analysis_names[s]]['simtype']=='md' 
					else 'meso_qs'):mscs[s].qmagst,
				('cgmd_t2d' if analysis_descriptors[analysis_names[s]]['simtype']=='md' 
					else 'meso_t2d'):mscs[s].t2d,
				'lenscale':msets[s].lenscale,
				}
		
#---comparison of curvature between MESO and CGMD methods
#---plots height-curvature correlation alongsize structure and variations

#---fixed limits for making smooth gifs
extremz_checkplot_global = -2.5,2.5
extrems_checkplot_global = 2.0
extremc0_checkplot_global = 0.03
#---calculations
mset = msets[0]
vecs = mean(mset.vecs,axis=0)
m,n = mset.griddims
#---figure
fig = plt.figure(figsize=(18,12))
gs = gridspec.GridSpec(5,2,hspace=0.15)
gs.update(left=0.0,right=0.3)
gs2 = gridspec.GridSpec(1,2)
gs2.update(left=0.4,right=1.0)
axeslist = []
#---plot curvature field
extrem = max([max([j.max() for j in mscs[i].c0s]) for i in range(2)])
if extremc0_checkplot_global != None: extrem = extremc0_checkplot_global
ax = plt.subplot(gs[0,0])
ax.imshow(mean(mscs[0].c0s,axis=0).T,vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
	interpolation='nearest',origin='lower')
axeslist.append(ax)
ax.set_title('CGMD')
ax = plt.subplot(gs[0,1])
axeslist.append(ax)
ax.set_title('MESO')
#---choose only a single frame, since the curvature field may be somewhat mobile at mesosscale
im = ax.imshow(mscs[1].c0s[0],vmax=extrem,vmin=0.,cmap=mpl.cm.binary,
	interpolation='nearest',origin='lower')
axins = inset_axes(ax,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,
	borderpad=0)
axeslist.append(axins)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
axins.set_ylabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',
	rotation=270)
#---plot average structure
vmax = max([mean(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
vmin = min([mean(msets[i].surf,axis=0).min()/msets[i].lenscale for i in range(2)])
extrem = max(abs(vmax),abs(vmin))
vmax,vmin = extrem,-extrem
if extremz_checkplot_global != None: vmin,vmax = extremz_checkplot_global
ax = plt.subplot(gs[1,0])
axeslist.append(ax)
im = ax.imshow(mean(msets[0].surf,axis=0).T/msets[0].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
	interpolation='nearest',origin='lower')
ax = plt.subplot(gs[1,1])
axeslist.append(ax)
im = ax.imshow(mean(msets[1].surf,axis=0).T/msets[1].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r,
	interpolation='nearest',origin='lower')
axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,borderpad=0)
axeslist.append(axins)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle (\mathrm{nm})$',rotation=270)
#---plot standard deviations
extrem = max([std(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
if extrems_checkplot_global != None: extrem = extrems_checkplot_global
ax = plt.subplot(gs[2,0])
axeslist.append(ax)
ax.imshow((std(msets[0].surf,axis=0)/msets[0].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.jet,
	interpolation='nearest',origin='lower')
ax = plt.subplot(gs[2,1])
axeslist.append(ax)
im = ax.imshow((std(msets[1].surf,axis=0)/msets[1].lenscale).T,vmin=0.,vmax=extrem,cmap=mpl.cm.jet,
	interpolation='nearest',origin='lower')
axins = inset_axes(ax,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,
	borderpad=0)
axeslist.append(axins)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
axins.set_ylabel(r'$\left\langle \left(z-\overline{z}\right)^{2} \right\rangle (\mathrm{{nm}^2})$',
	rotation=270)
#---compare 2D energy spectra
for m in [analysis_names.index(aname) for aname	in plot_reord]:
	mset = msets[m]
	axr2 = plt.subplot(gs[3,m])
	cm,cn = [int(i/2) for i in shape(mscs[m].tsum2d)]
	wid = 3
	dat = mscs[m].tsum2d[cm-wid:cm+wid+1,cn-wid:cn+wid+1]
	dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
	im = plotter2d(axr2,mset,dat=dat,
		cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
		ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
		fs=10,label_style='q',lims=[0.1,10],
		tickskip=int(round(mset.griddims[0]/6,-1)))
axins = inset_axes(axr2,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	borderpad=0)
axeslist.append(axins)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
axins.set_ylabel(
	r'$\left\langle \mathscr{H}_{el}\right\rangle \left(\frac{k_{B}T}{2}\right)^{-1}$',
	fontsize=fsaxlabel,rotation=270)
#---compare height-curvature correlation in 2D
for m in [analysis_names.index(aname) for aname	in plot_reord]:
	mset = msets[m]
	axr3 = plt.subplot(gs[4,m])
	cm,cn = [int(i/2) for i in shape(mscs[m].t2d[1])]
	wid = 3
	dat = mscs[m].t2d[1][cm-wid:cm+wid+1,cn-wid:cn+wid+1]
	dat[shape(dat)[0]/2,shape(dat)[1]/2] = (vmin+vmax)/2.
	im = plotter2d(axr3,mset,dat=dat,
		cmap=mpl.cm.RdBu_r,inset=False,cmap_washout=1.0,
		ticklabel_show=[1,1],tickshow=[1,1],centertick=False,
		fs=10,label_style='q',lims=[10**-5,10**-2],
		tickskip=int(round(mset.griddims[0]/6,-1)))
axins = inset_axes(axr3,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=axr3.transAxes,
	borderpad=0)
axeslist.append(axins)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
axins.set_ylabel(
	r'$\left\langle C_{0,q} h_{-q} \right\rangle $',
	fontsize=fsaxlabel,rotation=270)
axeslist.append(axins)
for ax in axeslist:
	plt.setp(ax.get_yticklabels(),fontsize=10)
	plt.setp(ax.get_xticklabels(),fontsize=10)	
#---spectrum plot
spectrum_summary(fig=fig,gs=gs2,
	titletext=r'$\mathrm{C_{0,hypo}}='+('{0:.3f}'.format(c0ask))+'a_0^{-1}'+\
	'='+('{0:.3f}'.format(c0ask*msets[1].lenscale))+'\:\mathrm{({nm}^{-1})}$')
#---save
plt.savefig(pickles+'fig-bilayer-couple-'+bigname+'.png',bbox_inches='tight',dpi=300)
if showplots: plt.show()
plt.close(fig)

'''
Note: make cool gifs with the following commands on light.site
gifmake ~/fig-bilayer-couple-compare-v614-120000-220000-200-v2005.gif 
	$(filefield fig-bilayer-couple-compare-v614-120000-220000-200-v2005-C_0-\*png C_0)
gifmake ~/fig-bilayer-couple-v614-120000-220000-200-v2005.gif 
	$(filefield fig-bilayer-couple-v614-120000-220000-200-v2005-C_0-\*png C_0)		
'''


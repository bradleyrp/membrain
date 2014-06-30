#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')
execfile('header-meso.py')

#---BATCH DEFAULTS
#-------------------------------------------------------------------------------------------------------------

batch_override = True

meso_avail = [
	'v2004',
	'v2005',
	'v2006',
	'v2009',
	'v2008',
	'v2012',
	'v2013',
	][-1:]

cgmd_avail = [
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	]
	
routine = [
	'calc',
	'masterplot',
	'checkplot',
	'plot2d',
	'plotphase'
	'print_new_c0_vals',
	][:2]

showplots = False
batch_override = True

p_interest = (meso_expt_toc[key])['parameter_name']

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def collapse_spectrum(qs,hs):
	'''For rectangular systems, treat lateral dimensions as equivalent and take their average.'''
	cen = array([i/2 for i in shape(qs)])
	discdists = [[sum(array([i-cen[0],j-cen[1]])**2) 
		for j in range(shape(qs)[1])] for i in range(shape(qs)[0])]
	uniqs = unique(flatten(discdists),return_inverse=True)[1]
	collapsed = [mean(ravel(hs)[where(uniqs==i)]) for i in range(uniqs.max()) 
		if mean(ravel(qs)[where(uniqs==i)]) != 0]
	return collapsed
	
#---LOAD SWEEP
#-------------------------------------------------------------------------------------------------------------

#---see if the coupling sweep has already been done
master_spectrum_dict = unpickle(
	pickles+'pkl.bilayer-coupling-sweep.'+'-'.join(meso_avail)+'-'+'-'.join(cgmd_avail)+'.pkl')

#---perform the sweep
if master_spectrum_dict == None:
	#---empty dictionary for cross-simulation comparisons
	master_spectrum_dict = {}
	#---this script will perform the script-coupling.py analysis for a parameter sweep
	for batch_cgmd in cgmd_avail:
		for batch_meso in meso_avail:
			#---set curvature extent from fixed parameter in the toc, assuming isotropic
			r_2 = (meso_expt_toc[batch_meso])['R_2']
			for c0ask in (meso_expt_toc[batch_meso])['parameter_sweep']: 
				match_scales = [batch_cgmd,batch_meso+'-'+p_interest+'-'+str(c0ask)]
				analysis_names = [batch_cgmd,batch_meso+'-'+p_interest+'-'+str(c0ask)]
				plot_reord = analysis_names
				bigname = '-'.join(analysis_names)
				if c0ask != 0:
					execfile('script-coupling.py')
					del msets,mscs,collect_c0s
	pickledump(master_spectrum_dict,
		'pkl.bilayer-coupling-sweep.'+'-'.join(meso_avail)+'-'+'-'.join(cgmd_avail)+'.pkl',
		directory=pickles)
else:
	routine = []
	execfile('script-coupling.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

def compute_undulation_residuals(c0asks,simnames,thresh=0.3,specnum=0,view_resids=False):
	'''Compute residuals between undulations for two simulations.'''
	subj = [[],[]]
	xdat = [[],[]]
	ydat = [[],[]]
	for s in range(len(simnames)):
		simname = simnames[s]
		#---use simulation indices to check CGMD vs MESO
		if simname != 'meso' and int(simname.split('-')[0].strip('v')) < 2000:
			subj[s] = master_spectrum_dict[simname+'-'+p_interest+'-'+str(c0asks[s])]
			xdat[s] = collapse_spectrum(subj[s]['cgmd_qs'],subj[s]['cgmd_qs'])
			ydat[s] = collapse_spectrum(subj[s]['cgmd_qs'],subj[s]['cgmd_t2d'][specnum])	
		#---perform mesoscale simulation lookup
		elif simname == 'meso':
			simname = [i for i in master_spectrum_dict.keys() if (float(i.split('-')[-1]) == c0asks[s] 
				and int(i.split('-')[0][1:]) > 2000)][0]
			subj[s] = master_spectrum_dict[simname]
			xdat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_qs'])
			ydat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_t2d'][specnum])
		else:
			raise Exception('except: no simulation found')
	dat0log = array([i for i in array([xdat[0],log10(ydat[0])]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat[1],log10(ydat[1])]).T if i[0] < thresh])
	resid = mean(abs(dat1log-dat0log)**2)
	tmp = dat1log[:,1]+mean(dat0log[:,1])-mean(dat1log[:,1])
	resid_shift = mean((dat0log-transpose([dat1log[:,0],tmp]))**2)
	if view_resids:
		ax = plt.subplot(111)
		plt.plot(dat0log[:,0],dat0log[:,1],'ro-')
		plt.plot(dat1log[:,0],dat1log[:,1],'bo-')
		plt.show()
	return resid,resid_shift
	
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
		#---plot lowest residual rankings
		tmp = array(comp_cgmd)[0].T
		for k in range(4):
			for i in range(len(tmp)):
				pos = argsort(tmp[i])[k]
				axr.text(i,pos,str(k),fontsize=16,horizontalalignment='center',va='center',weight='bold')
		#---settings
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
	plt.savefig(pickles+'fig-bilayer-couple-meta-'+filename_descriptor+\
		'-'.join([i.split('-')[0] for i in cgmd_avail]+meso_avail)+\
		'.png',
		bbox_inches='tight',dpi=300)
	if showplots: plt.show()
	plt.close(fig)

#---residual sweep parameters
cgmd_list = ['v614-120000-220000-200','v616-210000-310000-200']
qmagfilter = (analysis_descriptors[cgmd_avail[0]])['qmagfilter']
thresholds = [qmagfilter[1]]
thresholds = [0.2]
shift_curves = False

#---check the residuals to see that the wavevectors line up
if 0: compute_undulation_residuals([meso_c0s[0],meso_c0s[1]],['meso','meso'],
		thresh=qmagfilter[1],specnum=0,view_resids=True)
if 0: compute_undulation_residuals([meso_c0s[1],meso_c0s[1]],['meso',cgmd_avail[0]],
	thresh=qmagfilter[1],specnum=0,view_resids=True)

#---loop over comparisons
if 'comp' not in globals():
	#---perform the fit for both undulations and height-curvature correlations
	specnums = [0,1]
	fnames = ['undulations-','curvature_undulations-']	
	toptitles = ['Undulation residuals','Curvature-undulation residuals']
	#---loop
	comp = [[[] for i in thresholds] for sn in specnums]
	comp_cgmd = [[[] for i in thresholds] for sn in specnums]
	for sn in range(len(specnums)):	
		for thresh in thresholds:
			status('status: threshold = '+str(thresh))
			meso_c0s = list(sort(list(set([float(i.split('-')[-1]) 
				for i in master_spectrum_dict.keys() if int(i.split('-')[0][1:]) > 2000]))))
			residual_comparisons = zeros((len(meso_c0s),len(meso_c0s)))
			for i in meso_c0s:
				status('status: undulation comparison, c0 = '+str(i+1).ljust(6),
					i=meso_c0s.index(i),looplen=len(meso_c0s))
				for j in meso_c0s:
					resid,resid_shift = compute_undulation_residuals([i,j],['meso','meso'],
						thresh=thresh,specnum=sn)
					residual_comparisons[meso_c0s.index(i),meso_c0s.index(j)] = \
						(resid if not shift_curves else resid_shift)
			residual_comparisons_cgmd = zeros((len(meso_c0s),len(cgmd_list)))
			for i in meso_c0s:
				status('status: undulation comparison, c0 = '+str(i+1).ljust(4),
					i=meso_c0s.index(i),looplen=len(meso_c0s))
				for name in cgmd_list:
					resid,resid_shift = compute_undulation_residuals([i,i],['meso',name],
						thresh=thresh,specnum=sn)
					residual_comparisons_cgmd[meso_c0s.index(i),cgmd_list.index(name)] = \
						(resid if not shift_curves else resid_shift)
			comp[sn][thresholds.index(thresh)] = residual_comparisons
			comp_cgmd[sn][thresholds.index(thresh)] = residual_comparisons_cgmd
		#---get a0 from any mesoscale simulation
		simname = [i for i in master_spectrum_dict.keys() if (float(i.split('-')[-1]) == meso_c0s[0] 
			and int(i.split('-')[0][1:]) > 2000)][0]
		a0 = 1./(master_spectrum_dict[simname])['lenscale']
		plotter_undulation_residuals(thresholds,comp[sn],comp_cgmd[sn],toptitles[sn],a0,
			fnames[sn]+('shifted-' if shift_curves else ''))

#---generate a list of curvatures in the mesoscale simulations
if 'print_new_c0_vals' in routine:
	lister = [0.002, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02, 0.022, 
		0.025, 0.028, 0.03, 0.035, 0.04, 0.045, 0.05]
	lister = [0.002, 0.004, 0.006, 0.008, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 
		0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.032, 0.034, 
		0.036, 0.038, 0.04, 0.045, 0.05]
	print '( \''+'\' \''.join(['{0:.4f}'.format(round(i*a0,4)) for i in lister])+'\' )'



#!/usr/bin/python

'''
Notes:
	run this script after script-coupling-adv-batch.py
	if it can't find kappa_app then you need to rerun the sim
	
'''

from scipy.optimize import leastsq

qrange = [0.2,0.4]

def kappafits(scalefac):	
	subj = master_spectrum_dict[key]
	t2d = subj['cgmd_t2d']
	qmagst = subj['cgmd_qs']
	tsum2d = scalefac*(t2d[0]*qmagst**4-t2d[1]*qmagst**2-
		t2d[2]*qmagst**2+t2d[3])
	xdat = collapse_spectrum(qmagst,qmagst)
	ydat = collapse_spectrum(qmagst,tsum2d)
	#---target is global and shot is local, computed from the key each time
	#---? possibly very inefficient
	shot = array([xdat,ydat]).T
	if 0:
		#---screwed up!!!!!!!!!!!!!!!!!!!!!!!!!
		#---residual computation
		orders = [list(scipy.spatial.distance.cdist(
			[[i] for i in target.T[0]],[[i] for i in shot.T[0]]).argmin(axis=j)) for j in range(2)]
		inputs = [array([i for i in shot if i[0] < qrange[1] and i[0] > qrange[0]]),
			array([i for i in target if i[0] < qrange[1] and i[0] > qrange[0]])]
		if len(orders[0]) >= len(orders[1]): resid = sum((inputs[0].T[1]-inputs[1].T[1][orders[0]])**2)
		else: resid = sum((inputs[0].T[1][orders[1]]-inputs[1].T[1])**2)
		return resid
	inputs = [array([i for i in shot if i[0] < qrange[1] and i[0] > qrange[0]]),
		array([i for i in target if i[0] < qrange[1] and i[0] > qrange[0]])]
	return (mean(inputs[0].T[1])-mean(inputs[1].T[1]))**2

#---plotting the energy terms
if 1 or 'mesospec_collect' not in globals():

	mesospec_collect = []
	mesospec_collect_q = []
	meso_kappas = []
	cgmd_collect = []
	cgmd_collect_q = []
	cgmd_keys = []

	fig = plt.figure(figsize=(18,10))
	axes = []
	gs = gridspec.GridSpec(2,2)
	ax1 = fig.add_subplot(gs[0,0])
	axes.append(ax1)
	ax2 = fig.add_subplot(gs[0,1])
	axes.append(ax2)
	ax3 = fig.add_subplot(gs[1,0])
	axes.append(ax3)
	ax4 = fig.add_subplot(gs[1,1])
	axes.append(ax4)

	#---cgmd testnames for coloring
	cgmd_names = [analysis_descriptors[key]['testname'] 
		for key in analysis_descriptors.keys() if analysis_descriptors[key]['simtype']=='md']
	#---initial loop over all simulations loads them into tables for meso and CGMD
	for key in [i for i in master_spectrum_dict.keys() if i in analysis_descriptors.keys()]:
		status('status: plotting key = '+key)
		subj = master_spectrum_dict[key]
		
		#---rescale correlation terms to get the energy term
		scalefac = 1.
		if 'meso_qs' in subj.keys():
			t2d = subj['meso_t2d']
			qmagst = subj['meso_qs']
		elif 'cgmd_qs' in subj.keys():
			t2d = subj['cgmd_t2d']
			qmagst = subj['cgmd_qs']
		tsum2d = scalefac*(t2d[0]*qmagst**4-t2d[1]*qmagst**2-
			t2d[2]*qmagst**2+t2d[3])
		xdat = collapse_spectrum(qmagst,qmagst)
		ydat = collapse_spectrum(qmagst,tsum2d)
		
		#---save 1D spectra
		if 'meso_qs' in subj.keys(): mesospec_collect.append(ydat),mesospec_collect_q.append(xdat)
		if 'cgmd_qs' in subj.keys(): 
			cgmd_collect.append(ydat)
			cgmd_collect_q.append(xdat)
			cgmd_keys.append(key)
		
		#---plot specs and ordering by wavevector
		if (analysis_descriptors[key])['simtype'] == 'meso':
			color = 'k'
			alpha = 0.2
			zorder = 1
			lw = 1
			label = ''
			kappa_app = analysis_descriptors[key]['kappa_apparent']
		else:
			color = clrs[cgmd_names.index(analysis_descriptors[key]['testname'])]
			alpha = 0.0
			zorder = 2
			lw = 3
			label = analysis_descriptors[key]['label']
			mset = unpickle(pickles+analysis_descriptors[key]['locate'])
			kappa_app = mset.undulate_kappa
			print key+' kappa_apparent = '+str(kappa_app)
		meso_kappas.append(kappa_app)
		inds = argsort(array(xdat))
		
		#---CRITICAL PORTION OF THE CALCULATION
		'''
		Plan for each of the 1d kBT-equivalent spectra
			1. Choose a wavevector range to study
			2. Take the average of all mesoscale simulations, nearly equal (after kappa adjustment). 
			3. Then compute the supposed kappa for the CGMD simulations
			4. Then find the best fit to curvature
			5. Bask in the glory of a nice two-parameter fit?
		'''
		
		#---plot the raw spectra, withholding the CGMD here
		if (analysis_descriptors[key])['simtype'] == 'meso':
			ax1.plot(array(xdat)[inds],array(ydat)[inds]*kappa_app,
				color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)
			ax3.plot(array(xdat)[inds],array(ydat)[inds],
				color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)

		#---repeat the calculation with undulations only
		tsum2d = scalefac*(t2d[0]*qmagst**4)
		xdat = collapse_spectrum(qmagst,qmagst)
		ydat = collapse_spectrum(qmagst,tsum2d)
		#---plot undulations only
		if (analysis_descriptors[key])['simtype'] == 'meso':
			ax2.plot(array(xdat)[inds],array(ydat)[inds]*kappa_app,
				color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)
			ax4.plot(array(xdat)[inds],array(ydat)[inds],
				color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)

	#---plot the average mesoscale simulation after non-dimensionalizing by kappa	
	ax1.plot(mean(mesospec_collect_q,axis=0),mean(mesospec_collect,axis=0)*mean(meso_kappas),'k-',lw=4)
	ax3.plot(mean(mesospec_collect_q,axis=0),mean(mesospec_collect,axis=0),'k-',lw=4)
			
	#---rescale the CGMD simulations and determine kappa based on a wavevector magnitude range
	cgmd_keys = [key for key in analysis_descriptors.keys() if analysis_descriptors[key]['simtype'] == 'md']
	#---ACTUAL FITS
	for cnum in range(len(cgmd_keys)):
		key = cgmd_keys[cnum]
		bigref = array([mean(mesospec_collect_q,axis=0),mean(mesospec_collect,axis=0)*mean(meso_kappas)]).T
		target = bigref
		print 'optimizing'
		p_opt = leastsq(kappafits,1.)

		subj = master_spectrum_dict[key]
		scalefac = p_opt[0][0]
		print key+' scalefac = '+str(scalefac)
		t2d = subj['cgmd_t2d']
		qmagst = subj['cgmd_qs']
		tsum2d = scalefac*(t2d[0]*qmagst**4-t2d[1]*qmagst**2-
			t2d[2]*qmagst**2+t2d[3])
		xdat = collapse_spectrum(qmagst,qmagst)
		ydat = collapse_spectrum(qmagst,tsum2d)	
		color = clrs[cnum]
		alpha = 1.0
		zorder = 2
		lw = 3
		label = analysis_descriptors[key]['label']
		#---plot the rescaled CGMD simulation on the non-dimensionalized MESO curve
		ax1.plot(array(xdat)[inds],array(ydat)[inds],'-',
			color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)
		ax3.plot(array(xdat)[inds],array(ydat)[inds]/scalefac,'-',
			color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)

		tsum2d = t2d[0]*qmagst**4
		xdat = collapse_spectrum(qmagst,qmagst)
		ydat = collapse_spectrum(qmagst,tsum2d)
		ax2.plot(array(xdat)[inds],array(ydat)[inds]*scalefac,'-',
			color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)
		ax4.plot(array(xdat)[inds],array(ydat)[inds],'-',
			color=color,alpha=alpha,zorder=zorder,lw=lw,label=label)

	#---plot settings
	axes[0].set_title('energy term (rescaled by $\kappa$)')
	axes[1].set_title('undulations (rescaled by $\kappa$)')
	axes[2].set_title('energy term')
	axes[3].set_title('undulations')
	for ax in axes:
		ax.axvline(x=qrange[0],ymin=0,ymax=1.,lw=2,c='r',alpha=0.5,zorder=0)
		ax.axvline(x=qrange[1],ymin=0,ymax=1.,lw=2,c='r',alpha=0.5,zorder=0)
		ax.set_xscale('log')
		ax.set_yscale('log')
		h,l = ax.get_legend_handles_labels()
		h = [h[i] for i in range(len(l)) if l[i] != '']
		l = [l[i] for i in range(len(l)) if l[i] != '']
		ax.legend(h[::-1],l[::-1],loc='upper left')
	plt.savefig(pickles+'fig-bilayer-couple-meta-SPECTRA1.png',dpi=300)
	plt.show()


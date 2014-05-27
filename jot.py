#!/usr/bin/python

#---residuals BILAYER COUPLING

if 0:
	thresh = 1
	spec_query = [0]
	for s in range(len(spec_query)):
		#---original residuals
		a0 = msets[1].lenscale
		dat0 =array([i for i in mscs[0].t1d[spec_query[s]] if i[0] != 0 and i[0] < thresh])
		dat1 = array([i for i in mscs[1].t1d[spec_query[s]] if i[0] != 0 and i[0] < thresh])
		dat0log = log10(dat0)
		dat1log = log10(dat1)
		cd = scipy.spatial.distance.cdist(dat0log,dat1log)
		cd[cd.argmin(axis=argmax(shape(cd)))]
		resid = mean([cd[i,argmax(shape(cd))] for i in cd.argmin(axis=argmax(shape(cd)))])
		#---new residual signature curves
		datlogs = []
		for m in range(len(mscs)):
			xdat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
			datlogs.append([xdat0,ydat0])
		#---save and report
		print 'result: C_0 = '+('{0:.3f}'.format(c0ask*msets[1].lenscale)).rjust(5)+' (nm^-1)   '+\
			'resid = '+('{0:.3f}'.format(resid)).rjust(10)+'  status: npts = ('+str(len(dat0))+\
			','+str(len(dat1))+')'
		if 'collected_residuals' in globals(): 
			collected_residuals[cgmd_avail.index(cgmd_reference)][s].append([c0ask,resid])
			collected_residual_sigs[cgmd_avail.index(cgmd_reference)][s].append([c0ask,datlogs])


if 0:
	#---codeblock which compares undulations to their residuals, requires mscs
	ax1 = plt.subplot(211)
	ax1.plot(dat0[:,0],dat0[:,1],'r+')
	ax1.plot(dat1[:,0],dat1[:,1],'bx')
	m = 1
	xdat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
	ydat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
	ax.plot(xdat0,ydat0,'ro')
	m = 0
	xdat1 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
	ydat1 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
	ax1.plot(xdat1,ydat1,'bo')
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax2 = plt.subplot(212)
	dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
	resid_distn = cd.min(axis=argmax(shape(cd)))
	ax2.plot(dat0log[:,0],dat0log[:,1]-dat1log[:,1],'o-')
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xscale('log')
	plt.show()
	
if 0:	
	#---codeblock which looks at error signatures assuming same wavevectors
	ax2 = plt.subplot(111)
	dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
	ax2.plot(dat0log[:,0],dat0log[:,1]-dat1log[:,1],'o-')
	ax2.set_xscale('log')
	plt.show()
	
#---CONTINUING more advanced analysis of residuals BILAYER COUPLING

#---signatures of the residuals
if 0:
	thresh = 0.5
	fig = plt.figure(figsize=(12,12))
	ax2 = plt.subplot(111)
	for key in master_spectrum_dict.keys():
		subj = master_spectrum_dict[key]	
		xdat0 = collapse_spectrum(subj['cgmd_qs'],subj['cgmd_qs'])
		ydat0 = collapse_spectrum(subj['cgmd_qs'],subj['cgmd_t2d'][0])
		xdat1 = collapse_spectrum(subj['meso_qs'],subj['meso_qs'])
		ydat1 = collapse_spectrum(subj['meso_qs'],subj['meso_t2d'][0])
		dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
		dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
		resid = sum(abs(dat1log[:,1]-dat0log[:,1]))
		label = '{0:.2f}'.format(resid)
		ax2.plot(dat0log[:,0],dat1log[:,1]-dat0log[:,1],'-',label=label)
	ax2.legend(loc='upper right')
	ax2.set_xscale('log')
	ax2.set_xlim([min(dat0log[:,0])/2.,2.*max(dat0log[:,0])])
	plt.show()
	
	
if 0:
	pickledump(
		master_spectrum_dict,
		'pkl.bilayer-coupling-sweep.'+'-'.join(meso_avail)+'-'.join(cgmd_avail)+'.pkl',
		directory=pickles)
		
#---scatterplot of the mean residuals
if 0:
	#---prepare plot
	ax = plt.subplot(111)
	color_dict = {
		614:'r',
		616:'b',
		}
	thresholds = [0.3,0.9]
	alphas = [1.0,0.35]
	for thresh in thresholds:
		#---calculate
		residuals_cgmd_meso = []
		residuals_meso_meso = []
		cgmd_flags = []
		#subj_meso = master_spectrum_dict[('v614-120000-220000-200', 'v2005-C_0-0.038')]
		subj_meso = master_spectrum_dict[('v614-120000-220000-200', 'v2008-C_0-0.0147')]
		for key in master_spectrum_dict.keys():
			subj = master_spectrum_dict[key]	
			xdat0 = collapse_spectrum(subj['cgmd_qs'],subj['cgmd_qs'])
			ydat0 = collapse_spectrum(subj['cgmd_qs'],subj['cgmd_t2d'][0])
			xdat1 = collapse_spectrum(subj['meso_qs'],subj['meso_qs'])
			ydat1 = collapse_spectrum(subj['meso_qs'],subj['meso_t2d'][0])
			xdat2 = collapse_spectrum(subj_meso['meso_qs'],subj_meso['meso_qs'])
			ydat2 = collapse_spectrum(subj_meso['meso_qs'],subj_meso['meso_t2d'][0])
			dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
			dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
			dat2log = array([i for i in array([xdat2,log10(ydat2)]).T if i[0] < thresh])
			resid = mean(abs(dat1log[:,1]-dat0log[:,1]))
			resid_alt = mean(abs(dat1log[:,1]-dat2log[:,1]))
			label = '{0:.2f}'.format(resid)
			residuals_cgmd_meso.append([float(key[1].split('-')[-1])/2.32,resid])
			residuals_meso_meso.append([float(key[1].split('-')[-1])/2.32,resid_alt])
			flag = int(key[0].split('-')[0][1:])
			cgmd_flags.append(flag)
		#---plot
		for cgmd_flag in list(set(cgmd_flags)):
			plotdat = array([residuals_cgmd_meso[i] for i in range(len(residuals_cgmd_meso)) 
				if cgmd_flags[i] == cgmd_flag])
			ax.plot(plotdat[:,0],plotdat[:,1],'o',color=color_dict[cgmd_flag],
				alpha=alphas[thresholds.index(thresh)])
		cgmd_flag = 616
		plotdat = array([residuals_meso_meso[i] for i in range(len(residuals_meso_meso)) 
			if cgmd_flags[i] == cgmd_flag])
		ax.plot(plotdat[:,0],plotdat[:,1],'o',color='g',
			alpha=alphas[thresholds.index(thresh)])
	plt.show()


def 

if 'master_spectrum_dict' not in globals():
	master_spectrum_dict = unpickle(
		pickles+'pkl.bilayer-coupling-sweep.'+'-'.join(meso_avail)+'-'.join(cgmd_avail)+'.pkl')
	
	
	
	
	
	
	
	

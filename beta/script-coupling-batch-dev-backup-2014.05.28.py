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

#---plot the summary, see from script-coupling-batch.py now in jot.py for current development
if 0:
	annotate = False
	spec_colors = clrs[1:4]
	spec_labels = [r'$\left\langle h_{\mathbf{q}}h_{\mathbf{-q}}\right\rangle$',
		r'$\left\langle C_{0,\mathbf{q}} h_{-\mathbf{q}} \right\rangle $',
		r'$\left\langle C_{0,\mathbf{q}} C_{0,-\mathbf{q}} \right\rangle $']
	bigname_nometa = '-'.join(cgmd_avail+meso_avail)
	npanels = len(collected_residuals)
	fig = plt.figure(figsize=(10,2*npanels))
	gs = gridspec.GridSpec(npanels,5,hspace=0,wspace=0.1)
	for cri in range(npanels):
		ax = plt.subplot((gs[cri,:4] if npanels > 1 else gs[cri,:4])) 
		ax2 = plt.subplot((gs[cri,4] if npanels > 1 else gs[cri,4]),sharey = ax) 
		for spec_query in range(3):
			data = array([i for i in collected_residuals[cri][spec_query]])
			ax.scatter(data[:,0]/a0,data[:,1],s=40,color=spec_colors[spec_query],
				edgecolor='k',lw=1.)
			ax2.scatter(data[:,0]/a0,data[:,1],s=40,color=spec_colors[spec_query],
				edgecolor='k',lw=1.,
				label=(spec_labels[spec_query] if cri == 0 else None))
			if cri == 0 and spec_query == 2:
				ax2.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
		leftcut = 0.02
		ax.set_xlim(0,leftcut)
		label_offsets = [list(data[:,1]).index(i) for i in sort(data[:,1])]
		ax2.set_xlim(
			1/2.*(min([data[i,0]/a0 for i in label_offsets if data[i,0]/a0 >= leftcut])+\
			max([data[i,0]/a0 for i in label_offsets if data[i,0]/a0 < leftcut])),
			data[:,0].max()*1.1/a0)
		ax.spines['right'].set_visible(False)
		ax2.spines['left'].set_visible(False)
		ax.yaxis.tick_left()
		ax.tick_params(labeltop='off')
		ax2.yaxis.tick_right()
		plt.subplots_adjust(wspace=0.15)
		ax.set_yticklabels([])
		ax.set_ylim((0,array(collected_residuals)[...,1].max()*1.1))
		ax.grid(True)
		ax2.grid(True)
		ax.set_xticks(data[array([i for i in label_offsets if data[i,0]/a0 < leftcut]),0]/a0)
		ax2.set_xticks(data[array([i for i in label_offsets if data[i,0]/a0 >= leftcut]),0]/a0)
		if annotate:
			for ind in range(len(label_offsets)):
				x,y = data[label_offsets[ind]]
				plt.annotate(
					('{0:.3f}'.format(x)), 
					xy = (x/a0,y), 
					xytext = (100+(ind%2==1)*50,50+(ind%2==1)*20),
					textcoords = 'offset points', ha = 'right', va = 'bottom',
					bbox = dict(boxstyle='round,pad=0.2',fc = 'black', alpha = 0.2),
					arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
		if cri == 0:
			ax.set_title('CGMD-MESO matching, spectrum residuals')
			#ax2.legend(loc='center right',bbox_to_anchor=(1,0.5))
			#ax.legend(loc='upper right',bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)
		if cri == npanels-1:
			plt.setp(ax.xaxis.get_majorticklabels(),rotation=-90,fontsize=fsaxticks-4)
			plt.setp(ax2.xaxis.get_majorticklabels(),rotation=-90,fontsize=fsaxticks-4)
			ax.set_xlabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',fontsize=fsaxlabel)
		else:
			ax.set_xlabel(None)
			ax.set_xticklabels([])
			ax2.set_xlabel(None)
			ax2.set_xticklabels([])
		ax.set_ylabel((analysis_descriptors[cgmd_avail[cri]])['label'],fontsize=fsaxlabel)
	plt.savefig(pickles+'fig-bilayer-couple-meta-'+bigname_nometa+'.png',bbox_inches='tight')
	plt.show()

def compute_undulation_residuals(c0asks,simnames,thresh=0.3):
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
			ydat[s] = collapse_spectrum(subj[s]['cgmd_qs'],subj[s]['cgmd_t2d'][0])	
		#---perform mesoscale simulation lookup
		elif simname == 'meso':
			simname = [i for i in master_spectrum_dict.keys() 
				if (float(i.split('-')[-1]) == c0asks[s] and int(i.split('-')[0][1:]) > 2000)][0]
			subj[s] = master_spectrum_dict[simname]
			xdat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_qs'])
			ydat[s] = collapse_spectrum(subj[s]['meso_qs'],subj[s]['meso_t2d'][0])
		else:
			raise Exception('except: no simulation found')
	dat0log = array([i for i in array([xdat[0],log10(ydat[0])]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat[1],log10(ydat[1])]).T if i[0] < thresh])
	resid = mean(abs(dat1log-dat0log)**2)
	return resid

#---residual sweep parameters
cgmd_list = ['v614-120000-220000-200','v616-210000-310000-200']
thresholds = [0.3]

if 'comp' not in globals():
	comp = [[] for i in thresholds]
	comp_cgmd = [[] for i in thresholds]  
	for thresh in thresholds:
		status('status: threshold = '+str(thresh))
		meso_c0s = list(sort(list(set([float(i.split('-')[-1]) 
			for i in master_spectrum_dict.keys() if int(i.split('-')[0][1:]) > 2000]))))
		residual_comparisons = zeros((len(meso_c0s),len(meso_c0s)))
		for i in meso_c0s:
			status('status: undulation comparison, c0 = '+str(i+1).ljust(6),i=meso_c0s.index(i),looplen=len(meso_c0s))
			for j in meso_c0s:
				residual_comparisons[meso_c0s.index(i),meso_c0s.index(j)] = \
					compute_undulation_residuals([i,j],['meso','meso'],thresh=thresh)
		residual_comparisons_cgmd = zeros((len(meso_c0s),len(cgmd_list)))
		for i in meso_c0s:
			status('status: undulation comparison, c0 = '+str(i+1).ljust(4),i=meso_c0s.index(i),looplen=len(meso_c0s))
			for name in cgmd_list:
				residual_comparisons_cgmd[meso_c0s.index(i),cgmd_list.index(name)] = \
					compute_undulation_residuals([i,i],['meso',name],thresh=thresh)
		comp[thresholds.index(thresh)] = residual_comparisons
		comp_cgmd[thresholds.index(thresh)] = residual_comparisons_cgmd

fig = plt.figure(figsize=((4,3*len(thresholds)) if len(thresholds)>1 else (8,8)))
gs = gridspec.GridSpec(len(thresholds),2,wspace=0.1,hspace=0.1,
	width_ratios=[len(meso_c0s),len(cgmd_list)])
cmap = mpl.cm.jet
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
	axl.set_xticklabels(['{:.3f}'.format(i) for i in meso_c0s])
	axl.set_yticklabels(['{:.3f}'.format(i) for i in meso_c0s])
	for label in im.axes.xaxis.get_ticklabels():
		label.set_rotation(90)
	axl.set_xlim((0,len(meso_c0s)))
	im = axr.imshow(residual_comparisons_cgmd,interpolation='nearest',cmap=cmap,
		vmax=max_resid,vmin=0)
	axr.set_xticks(range(len(cgmd_list)))
	axr.set_xticklabels([(analysis_descriptors[name])['label'] for name in cgmd_list])
	for label in im.axes.xaxis.get_ticklabels():
		label.set_rotation(90)
	plt.setp(axr.get_yticklabels(), visible=False)
	axins = inset_axes(axr,width=1./len(cgmd_list)/2.,height="100%",loc=3,
		bbox_to_anchor=(1.2,0.0,1.,1.),
		bbox_transform=axr.transAxes,
		borderpad=0)
	ticks = list(linspace(0,max_resid,21))
	norm = mpl.colors.BoundaryNorm(ticks,cmap.N)
	cbar = mpl.colorbar.ColorbarBase(axins,cmap=cmap,ticks=ticks,boundaries=ticks,
		norm=norm,spacing='proportional',format='%03f')
	cbar.set_ticklabels([str(i) for i in ticks])
	axl.set_xlim((-0.5,len(meso_c0s)-1+0.5))
	axr.set_xlim((-0.5,len(cgmd_list)-1+0.5))
	axl.set_title('Undulation Residuals',fontsize=fsaxtitle)
	axl.set_ylabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
	axl.set_xlabel(r'$\mathrm{mesoscale\:C_0\:({nm}^{-1})}$',fontsize=fsaxlabel)
plt.show()
	

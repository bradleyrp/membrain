#!/usr/bin/python

#def view_figures_timeseries(params=None,maxhs=None,maxhxys=None):
if 1:
	#---Calculate point-wise residuals
	#resids = [mean([abs(gauss2d(params[0],i[0],i[1])-i[2]) for i in target_zones[j]]) for j in range(len(params))]
	#nfitpts = [len(i) for i in target_zones]
	mpl.rcParams.update({'font.size': 16})
	fig = plt.figure(figsize=(10,10))
	ax1 = plt.subplot2grid((4,1),(0,0))
	ax1.plot(range(len(validhis)),validhs)
	ax1.set_xlim((0,len(validhis)))
	ax1.set_ylabel('$H_{max}$')
	ax1.grid()
	plt.title('Framewise measurements and residuals')
	ax2 = plt.subplot2grid((4,1),(1,0))
	ax2.plot(range(len(validhis)),sigma_x)
	ax2.plot(range(len(validhis)),sigma_y)
	ax2.set_yscale('log')
	ax2.set_xlim((0,len(validhis)))
	ax2.set_ylabel('($\sigma_x,\sigma_y)$')
	ax2.grid()
	ax3 = plt.subplot2grid((4,1),(2,0))
	ax3.plot(range(len(validhis)),[resids[i] for i in validhis])
	ax3.set_xlim((0,len(validhis)))
	ax3.set_ylabel('$\sum(z-\hat{z})$')
	ax3.grid()
	ax4 = plt.subplot2grid((4,1),(3,0))
	ax4.plot(range(len(validhis)),[nfitpts[i] for i in validhis])
	ax4.set_xlim((0,len(validhis)))
	ax4.set_ylabel('$N_{points}$')
	ax4.set_xlabel('Frame')
	ax4.grid()
	#plt.savefig(pickles+'result.fig.dimple.'+sigma_ytemprefix+'framewise'+'.png',dpi=300,bbox_extra_artists=[ylabel1,ylabel2],bbox_inches='tight')
	#plt.close()
	plt.show()

#def view_figures_curvatures(params=None,maxhs=None,maxhxys=None,scales='linear',fullrange=False,extrasuffix=''):
if 0:
	scales='linear'
	fullrange = False
	extrasuffix = ''

	if fullrange == True:
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
	else:
		validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] > curvature_filter[0] 
			and 10*maxhs[i] < curvature_filter[1])]
	#---Nanometer correction in the calculation
	validhs = [10*maxhs[i] for i in validhis]
	mpl.rcParams.update({'font.size': 14})
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((2,1),(0,0))
	#---Nanometer correction in the calculation
	if scales == 'log-transformed':
		ax.hist([log10(i) for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8)
	elif scales == 'linear':
		ax.hist([i for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,bins=40)		
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	if scales == 'log-transformed':
		mean_validhs = exp(mean(log(validhs)))
		plt.xlabel('Log $H_{max}$ $(nm^{-1})$', labelpad = 10)
	elif scales == 'linear':
		mean_validhs = mean(validhs)
		plt.xlabel(r'$H_{max}\:(\textsf{nm}^{-1})$', labelpad = 10)	
	ax.text(0.65,0.75,r'$\left\langle H_{max}\right\rangle='+str('%.3f'%mean_validhs)+
		r'\:\textsf{nm}^{-1}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+
		'$ frames',transform=ax.transAxes,fontsize=14)
	ylabel1 =plt.ylabel('Frames', labelpad = 10)
	plt.title('Fitted $H_{max}$, within $'+str(cutoff_distance)+'$ nm of protein')
	legend1 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	ax.grid()
	ax2 = plt.subplot2grid((2,1),(1,0))
	#---Nanometer correction below
	if scales == 'log-transformed':
		sigma_x = [log10(abs(params[i][4])/10.) for i in validhis]
		sigma_y = [log10(abs(params[i][5])/10.) for i in validhis]
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.6,histtype='stepfilled')
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.6,histtype='stepfilled')
		
	elif scales == 'linear':
		sigma_x = [abs(params[i][4])/10. for i in validhis]
		sigma_y = [abs(params[i][5])/10. for i in validhis]
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.6,histtype='stepfilled',bins=250)
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.6,histtype='stepfilled',bins=250)
		ax2.set_xlim((0,50))
		print 'Extent data are really skewed. Consider norming this one.'
		print 10**mean([log10(abs(params[i][4])/10.) for i in validhis])
		print 10**mean([log10(abs(params[i][5])/10.) for i in validhis])
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.spines['bottom'].set_position(('outward', 20))
	ax2.spines['left'].set_position(('outward', 30))
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	if scales == 'log-transformed':
		plt.xlabel('Log extent $(\sigma_a,\sigma_b)$ $\log_{10}$(nm)', labelpad = 10)
		ax2.text(0.05,0.65,r'$\sigma_b='+str('%.3f'%exp(mean(sigma_y)))+'$ $\\textsf{nm}$',transform=ax2.transAxes,fontsize=14)
	elif scales == 'linear':
		plt.xlabel(r'Extent $(\sigma_a,\sigma_b)\:(\textsf{nm})$', labelpad = 10)
		ax2.text(0.65,0.75,r'$\sigma_a='+str('%.3f'%mean(sigma_x))+r'\:\textsf{nm}$'+'\n'+r'$\sigma_b='+
			str('%.3f'%mean(sigma_y))+r'\:\textsf{nm}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+
			'$ frames',transform=ax2.transAxes,fontsize=14)
	plt.title('Fitted extent ($\sigma_a,\sigma_b$), within $'+str(cutoff_distance)+'$ nm of protein')
	ylabel2 = plt.ylabel('Frames', labelpad = 10)
	legend2 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	ax2.grid()
	plt.tight_layout()
	#plt.savefig(pickles+'result.fig.dimple.'+sigma_ytemprefix+extrasuffix+'.png',dpi=300,bbox_extra_artists=[ylabel1,ylabel2],bbox_inches='tight')
	#plt.close()
	plt.show()
	


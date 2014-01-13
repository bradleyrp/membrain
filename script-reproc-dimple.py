#!/usr/bin/python -i

if 1:

	from membrainrunner import *

	import os
	from scipy.optimize import leastsq
	import matplotlib as mpl
	from pylab import *
	#mpl.rc('text.latex', preamble='\usepackage{sfmath}')
	#mpl.rcParams['axes.linewidth'] = 2.0

	skip = 1
	framecount = None
	location = ''
	execfile('locations.py')
	
	#---analysis plan
	analysis_plan = slice(-2,-1)
	analysis_descriptors = [
		('pkl.dimple.v614-stress.md.part0002.rerun.pkl'),
		()]

		
	
	result_data_collection = pickle.load(open(pickles+'pkl.dimple.v614-stress.md.part0002.rerun.pkl','r'))

	which_brewer_colors = [0,1]
	clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]

if 1:
	means='linear'
	scales='linear'
	fullrange=True
	extrasuffix=''
	sigma_calc='mode'


	'''
	possible plots
	
	1. all ENTH systems on one plot
	2. all ENTH systems on a plot stacked above Exo70 systems
	3. all systems on one plot
	
	'''
	single_plot = 1
	if single_plot:

		params_neg = result_data_collection[0].get(['type','params'])
		maxhs_neg = result_data_collection[0].get(['type','maxhs'])
		maxhxys_neg = result_data_collection[0].get(['type','maxhxys'])

		params_poz = result_data_collection[1].get(['type','params'])
		maxhs_poz = result_data_collection[1].get(['type','maxhs'])
		maxhxys_poz = result_data_collection[1].get(['type','maxhxys'])
	
		params_batch = (params_neg,params_poz)
		maxhs_batch = (maxhs_neg,maxhs_poz)
		maxhxys_batch = (maxhxys_neg,maxhxys_poz)

		'''
		params = (params_neg,params_poz)
		maxhs = (maxhs_neg,maxhs_poz)
		maxhxys = (maxhxys_neg,maxhxys_poz)

		(params_neg,params_poz) = params
		(maxhs_neg,maxhs_poz) = maxhs
		(maxhxys_neg,maxhxys_poz) = maxhxys
		'''
	
		#---top plot
		fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(10,4))
		nbins = 20
		minval,maxval = -0.10,0.10
		
		params, maxhs, maxhxys = params_neg, maxhs_neg, maxhxys_neg
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
		#---nanometer correction
		validhs = [10*maxhs[i] for i in validhis]
		hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		axes.plot(mid0,hist0,c=clrs[1],alpha=1.,lw=2)
		axes.fill_between(mid0,hist0,[0 for i in mid0],facecolor=clrs[1],alpha=0.35,interpolate=True)
		
		params, maxhs, maxhxys = params_poz, maxhs_poz, maxhxys_poz
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
		#---nanometer correction
		validhs = [10*maxhs[i] for i in validhis]
		#---plot
		hist0,binedge0 = numpy.histogram(validhs,bins=nbins,normed=True,range=(minval,maxval))
		mid0 = (binedge0[1:]+binedge0[:-1])/2
		axes.plot(mid0,hist0,c=clrs[0],alpha=1.,lw=2)
		axes.fill_between(mid0,hist0,[0 for i in mid0],facecolor=clrs[0],alpha=0.35,interpolate=True)
	


		plt.show()

if 0:
	if fullrange == True:
		#---Nanometer correction
		validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > 10**-5 and abs(10*maxhs[i]) < 0.1)]
	else:
		#---Nanometer correction
		if height_direction == 0:
			print 'hdir 0'
			validhis = [i for i in range(len(maxhs)) if (abs(10*maxhs[i]) > curvature_filter[0] 
				and abs(10*maxhs[i]) < curvature_filter[1])]
		elif height_direction == 1:
			print 'hdir 1'
			validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] > curvature_filter[0] 
				and 10*maxhs[i] < curvature_filter[1])]
		elif height_direction == -1:
			print 'hdir -1'
			validhis = [i for i in range(len(maxhs)) if (10*maxhs[i] < -1*curvature_filter[0] 
				and 10*maxhs[i] > -1*curvature_filter[1])]
	#---Nanometer correction
	validhs = [10*maxhs[i] for i in validhis]
	#mpl.rcParams.update({'font.size': 30})
	#mpl.rcParams.update({'font.style':'sans-serif'})
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((2,1),(0,0))
	#---Nanometer correction in the calculation
	if scales == 'log':
		ax.hist([log10(i) for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,linewidth=2)
	elif scales == 'linear':
		ax.hist([i for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,bins=40,linewidth=2)		
	#ax.spines['top'].set_visible(False)
	#ax.spines['right'].set_visible(False)
	#ax.spines['bottom'].set_position(('outward', 20))
	#ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	#print ax.get_ticks()
	#ax.yaxis.set_ticks([i for i in ax.yaxis.get_ticks() if (i%100) == 0])
	ax.set_xlim((-0.10,0.10))
	if means == 'linear':
		mean_validhs = mean(validhs)
	elif means == 'log':
		mean_validhs = exp(mean(log(validhs)))
	if scales == 'log':
		plt.xlabel('Log $H_{max}$ $(nm^{-1})$', labelpad = 10)
	elif scales == 'linear':
		plt.xlabel(r'$H_{max}\:(\textsf{nm}^{-1})$', labelpad = 10)	
	#ax.text(0.65,0.75,r'$\left\langle H_{max}\right\rangle='+str('%.3f'%mean_validhs)+
	#	r'\:\textsf{nm}^{-1}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+'$ frames',
	#	transform=ax.transAxes,fontsize=28)
	ylabel1 =plt.ylabel('Frames', labelpad = 10)
	#plt.title('Fitted $H_{max}$, within $'+str(cutoff_distance)+'$ nm of protein')
	legend1 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax.grid()
	ax2 = plt.subplot2grid((2,1),(1,0))
	#---Nanometer correction below
	if means == 'linear':
		sigma_x = [abs(params[i][4])/10. for i in validhis]
		sigma_y = [abs(params[i][5])/10. for i in validhis]
	elif means == 'log':
		sigma_x = [log10(abs(params[i][4])/10.) for i in validhis]
		sigma_y = [log10(abs(params[i][5])/10.) for i in validhis]
	if scales == 'log':
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.8,histtype='stepfilled',linewidth=2)
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.8,histtype='stepfilled',linewidth=2)
	elif scales == 'linear':
		histplot2 = ax2.hist(sigma_x,color=colorcodes[1],alpha=0.8,histtype='stepfilled',bins=250,linewidth=2)
		histplot1 = ax2.hist(sigma_y,color=colorcodes[2],alpha=0.8,histtype='stepfilled',bins=250,linewidth=2)
		ax2.set_xlim((0,50))
		print 'Extent data are really skewed. Consider norming this one.'
		print 10**mean([log10(abs(params[i][4])/10.) for i in validhis])
		print 10**mean([log10(abs(params[i][5])/10.) for i in validhis])
	#ax2.spines['top'].set_visible(False)
	#ax2.spines['right'].set_visible(False)
	#ax2.spines['bottom'].set_position(('outward', 20))
	#ax2.spines['left'].set_position(('outward', 30))
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	for tick in ax2.xaxis.get_major_ticks():
		tick.label.set_fontsize(22)
	for tick in ax2.yaxis.get_major_ticks():
		tick.label.set_fontsize(22) 
	if sigma_calc == 'mean':
		sigma_a_calc = mean(sigma_x)
		sigma_b_calc = mean(sigma_y)
	elif sigma_calc == 'logmean':
		sigma_a_calc = exp(mean(sigma_x))
		sigma_b_calc = exp(mean(sigma_y))
	elif sigma_calc == 'mode':
		hist1 = histogram(sigma_x,bins=250)
		hist2 = histogram(sigma_y,bins=250)
		sigma_a_calc = hist1[1][argmax(hist1[0])]
		sigma_b_calc = hist2[1][argmax(hist2[0])]
	if scales == 'log-transformed':
		plt.xlabel('Log extent $(\sigma_a,\sigma_b)$ $\log_{10}$(nm)', labelpad = 10)
		#ax2.text(0.05,0.65,r'$\sigma_a='+str('%.3f'%sigma_a_calc)+'$ $\\textsf{nm}$'+'\n'+r'$\sigma_b='+
		#	str('%.3f'%sigma_b_calc)+'$ $\\textsf{nm}$',transform=ax2.transAxes,fontsize=14)
	elif scales == 'linear':
		plt.xlabel(r'Extent $(\sigma_a,\sigma_b)\:(\textsf{nm})$', labelpad = 10)
		#ax2.text(0.65,0.75,r'$\sigma_a='+str('%.3f'%sigma_a_calc)+r'\:\textsf{nm}$'+'\n'+r'$\sigma_b='+
		#	str('%.3f'%sigma_b_calc)+r'\:\textsf{nm}$'+'\nfitted $'+str(len(validhs))+'/'+str(len(maxhs))+
		#	'$ frames',transform=ax2.transAxes,fontsize=14)
	#plt.title('Fitted extent ($\sigma_a,\sigma_b$), within $'+str(cutoff_distance)+'$ nm of protein')
	ylabel2 = plt.ylabel('Frames', labelpad = 10)
	legend2 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={'size':22})
	ax2.grid()
	plt.tight_layout()
	plt.savefig(pickles+'fig-'+sysname+'-dimple.png',dpi=500,
		bbox_extra_artists=[ylabel1,ylabel2],bbox_inches='tight')
	plt.close()


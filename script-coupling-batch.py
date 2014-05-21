#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')
execfile('header-meso.py')

meso_avail = [
	'v2004',
	'v2005',
	'v2006',
	]
cgmd_avail = [
	'v614-120000-220000-200',
	'v616-210000-310000-200',
	]

#---prepare list of collected residuals for a sweep over available simulations
#collected_residuals = [[[],[],[]] for i in range(len(cgmd_avail))]
#collected_residuals_sigs = [[[],[],[]] for i in range(len(cgmd_avail))]
#thresh = 0.3

#---empty dictionary for cross-simulation comparisons
master_spectrum_dict = {}

#---this script will perform the script-coupling.py analysis for a parameter sweep
for batch_cgmd in cgmd_avail:
	for batch_meso in meso_avail:
		for c0ask in (meso_expt_toc[batch_meso])['parameter_sweep']: 
			execfile('script-coupling.py')
			del msets,mscs,collect_c0s

#---plot the summary
if collected_residuals != [] and 0:
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

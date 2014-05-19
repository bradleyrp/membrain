#!/usr/bin/python

if 0:
	from mayavi import mlab
	meshplot(collect_c0s[0][0],show='surf')
	meshplot(collect_c0s[1][0],show='surf')
if 0:
	print array(mscs[0].c0s).max()
	print array(mscs[1].c0s).max()

	maxim = max([array(mscs[0].c0s).max(),array(mscs[1].c0s).max()])

	ax = plt.subplot(121)
	ax.imshow(mscs[0].c0s[0],vmin=0,vmax=maxim)
	ax = plt.subplot(122)
	ax.imshow(mscs[1].c0s[0],vmin=0,vmax=maxim)
	plt.show()

if 0:
#---plot the summary
#if collected_residuals != []:
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
		leftcut = 0.1
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
	
'''
>>> a0
0.40690525690439677
>>> array([0.008,0.01,0.015,0.02,0.025,0.030,0.035,0.04,0.045,0.05])/a0
array([ 0.0196606 ,  0.02457575,  0.03686362,  0.04915149,  0.06143936,
        0.07372724,  0.08601511,  0.09830298,  0.11059085,  0.12287873])
'''

def collapse_spectrum_prev(subj):
	cen = array([i/2 for i in shape(subj)])
	return [mean([
		subj[cen[0]][cen[1]+cy],
		subj[cen[0]][cen[1]-cy],
		subj[cen[0]+cx][cen[1]],
		subj[cen[0]-cx][cen[1]],
		]) for cy in range(1,cen[1])
		for cx in range(1,cen[0])]


#---debug
if 0:
	subj = mscs[1].qmagst
	cen = array([i/2 for i in shape(subj)])
	answq = [[
		subj[cen[0]][cen[1]+cy],
		subj[cen[0]][cen[1]-cy],
		subj[cen[0]+cx][cen[1]],
		subj[cen[0]-cx][cen[1]],
		] for cx in range(1,cen[0])
		for cy in range(1,cen[1])]
	subj = mscs[1].t2d[0]
	answh = [[
		subj[cen[0]][cen[1]+cy],
		subj[cen[0]][cen[1]-cy],
		subj[cen[0]+cx][cen[1]],
		subj[cen[0]-cx][cen[1]],
		] for cx in range(1,cen[0])
		for cy in range(1,cen[1])]
if 0:
	print answq[1]
	print answh[1]
	#print collapse_spectrum(mscs[1].qmagst)[0]
	#print collapse_spectrum(mscs[1].t2d[0])[0]
	
#---reformulate

def collapse_spectrum(subj):
	cen = array([i/2 for i in shape(subj)])
	return [mean([
		subj[cen[0]][cen[1]+cy],
		subj[cen[0]][cen[1]-cy],
		subj[cen[0]+cx][cen[1]],
		subj[cen[0]-cx][cen[1]],
		]) for cy in range(1,cen[1])
		for cx in range(1,cen[0])]
		
if 1:
	discgrid = [[cen[0]+cy+i[0],cen[1]+cx+i[1]] for i in [[0,1],[0,-1],[1,0],[-1,0]] for cy in range(1,cen[1]) for cx in range(1,cen[0])]
	discgrid = [[sqrt(sum(array([cx,cy])**2)) for cy in range(-cen[1]+1,cen[1]+1)] for cx in range(-cen[0]+1,cen[0]+1)]
	plt.imshow(array(discgrid).T,interpolation='nearest',origin='lower')
	plt.show()
	tmp = unique(discgrid,return_inverse=True)[1]


#---demonstrate
if 0:
	qtmp = mscs[1].t1d[0][:,0]
	htmp = mscs[1].t1d[0][:,1]
	tmp = unique(qtmp,return_inverse=True)[1]
	htmp_reduce = [mean(htmp[where(tmp==i)]) for i in range(tmp.max())]
	qtmp_reduce = [mean(qtmp[where(tmp==i)]) for i in range(tmp.max())]
	ax = plt.subplot(111)
	ax.plot(qtmp,htmp,'bo')
	if 0: ax.plot(qtmp_reduce,htmp_reduce,'rx')
	ax.plot(collapse_spectrum(mscs[1].qmagst),collapse_spectrum(mscs[1].t2d[0]),'m+')
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.show()





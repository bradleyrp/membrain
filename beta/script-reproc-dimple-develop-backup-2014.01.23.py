#!/usr/bin/python -i

v1 = False
v2 = True

#---development attempt 2 here, where I write code to plot visual summaries of the fits from results_stack
if v2:

	from scipy import ndimage

	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

	analysis_descriptors = [('pkl.dimple.v614-stress.md.part0002.rerun.pkl',(clrs[0],clrs[1]),
		r'$\textbf{{ENTH}\ensuremath{\times}4}$',1)]

	ad = analysis_descriptors[0]
	resdat = pickle.load(open(pickles+subdir+ad[0],'r'))
	mset = unpickle(pickles+resdat[0].getnote('startpickle'))	
	
	fr = 10
	vecs = mset.vec(fr)
	griddims = mset.griddims
	cutoff_distance = resdat[0].getnote('cutoff_distance')
	cutoff = cutoff_distance*10/(vecs[0]/griddims[0])
	zfilterdir = 1
	show_colorbar = False
	
	#---tesselate (via PBCs) the original surface
	#---surfpbc = numpy.tile(mset.surf[fr],(3,3))
	surfpbc = mset.surf[fr]
	#---positive/negative domain selection
	surfpoz = array([[(1 if surfpbc[i][j] > 0 else 0) 
		for j in range(1*(griddims[1]-1))] for i in range(1*(griddims[0]-1))]).T
	surfneg = array([[(1 if surfpbc[i][j] < 0 else 0) 
		for j in range(1*(griddims[1]-1))] for i in range(1*(griddims[0]-1))]).T
	label_im_poz, nb_labels = ndimage.label(surfpoz)
	label_im_neg, nb_labels = ndimage.label(surfneg)
	surf_discrete = label_im_poz - label_im_neg
	#---protein shadow selection
	prot_disc_pts = [(int(round(pt[0]/vecs[0]*(griddims[0]-1))),int(round(pt[1]/vecs[1]*(griddims[1]-1)))) 
		for pt in mset.protein[fr]]
	#---the following is equivalent
	#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	prot_locs = [[(1 if (i,j) in prot_disc_pts else 0) for i in range(griddims[0]-1)] 
		for j in range(griddims[1]-1)]
	#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	gridpos = [[i,j] for i in range(griddims[0]-1) for j in range(griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(prot_disc_pts,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	shadow = [[(1 if [i,j] in bufferlist else 0) for i in range(griddims[0]-1)] 
		for j in range(griddims[1]-1)]
	#---final selection
	target = [list(i) for i in array(where(((array(surf_discrete).T>0.)+(array(shadow).T==1)))).T]
	target = [list(i) for i in array(where((array(surf_discrete).T>0.)*(array(shadow).T==1)!=0)).T]
	if zfilterdir == 1:
		supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in target])) if i > 0]
	else:
		supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in target])) if i < 0]
	fintarget = [list([l[1],l[0]]) for l in [j for k in [list(array(where(surf_discrete==i)).T) 
		for i in supplement_domains] for j in k]]
	#---convert to xyz
	targetxyz = array([[i[0]*vecs[0]/(griddims[0]-1),i[1]*vecs[1]/(griddims[1]-1),
		mset.surf[fr][i[0],i[1]]] for i in fintarget])
		
	#---calculate fitted height map
	zmap = array([[gauss2d(params[fr],x*vecs[0]/(griddims[0]-1),y*vecs[1]/(griddims[1]-1)) 
		for y in range(griddims[1]-1)] for x in range(griddims[0]-1)]).T
	maxzpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
		for j in [targetxyz[argmax([gauss2d(params[fr],i[0],i[1]) for i in targetxyz])]]][0]
	minzpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
		for j in [targetxyz[argmin([gauss2d(params[fr],i[0],i[1]) for i in targetxyz])]]][0]
	minz = np.min([0.1*gauss2d(params[fr],i[0],i[1]) for i in targetxyz])
	maxz = np.max([0.1*gauss2d(params[fr],i[0],i[1]) for i in targetxyz])

	#---calculate curvature map
	curvmap = array([[-1*gauss2dh(params[fr],x*vecs[0]/(griddims[0]-1),y*vecs[1]/(griddims[1]-1)) 
		for y in range(griddims[1]-1)] for x in range(griddims[0]-1)]).T
	maxhpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
		for j in [targetxyz[argmax([-1*gauss2dh(params[fr],i[0],i[1]) for i in targetxyz])]]][0]
	minhpos = [(j[0]/vecs[0]*(griddims[0]-1),j[1]/vecs[1]*(griddims[1]-1)) 
		for j in [targetxyz[argmin([-1*gauss2dh(params[fr],i[0],i[1]) for i in targetxyz])]]][0]
	minh = np.min([-10*gauss2dh(params[fr],i[0],i[1]) for i in targetxyz])
	maxh = np.max([-10*gauss2dh(params[fr],i[0],i[1]) for i in targetxyz])

	#---calculate residual map
	residmap = zmap-mset.surf[fr][:-1,:-1]
	maxresid = 0.1*np.max(residmap)
	minresid = 0.1*np.min(residmap)

	fig = plt.figure(figsize=(8,8*0.965))
	gs = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0)
	domaincountlim = max(np.max(surf_discrete),abs(np.min(surf_discrete)))
	ax = fig.add_subplot(gs[0,0])
	img = ax.imshow(surf_discrete,extent=None,origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-domaincountlim,vmax=domaincountlim)
	for pt in fintarget:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
	for pt in bufferlist:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./5/100*griddims[0],color='w',alpha=1.))
	for pt in [list(i) for i in array(where(array(prot_locs)==1)).T]:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./3*2/100*griddims[0],color='k',alpha=1.))
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_title('DOMAINS',fontsize=16)

	ax = fig.add_subplot(gs[1,0])
	curvmax = np.max(abs(curvmap))
	img = ax.imshow(curvmap,extent=None,origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-curvmax,vmax=curvmax)
	for pt in fintarget:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
	ax.add_patch(plt.Circle(maxhpos,radius=2./100*griddims[0],color='k',alpha=1.))
	ax.add_patch(plt.Circle(minhpos,radius=2./100*griddims[0],color='k',alpha=1.))
	ax.add_patch(plt.Circle(maxhpos,radius=1./100*griddims[0],color='r',alpha=1.))
	ax.add_patch(plt.Circle(minhpos,radius=1./100*griddims[0],color='b',alpha=1.))
	if show_colorbar:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img,cax=cax)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_xlabel(r'$\mathrm{H\,({nm}^{-1})\in['+str('%3.3f'%minh)+','+str('%3.3f'%maxh)+']}$',fontsize=16)

	ax = fig.add_subplot(gs[1,1])
	zmax = np.max(abs(zmap))
	img = ax.imshow(zmap,extent=None,origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-zmax,vmax=zmax)
	for pt in fintarget:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
	ax.add_patch(plt.Circle(maxzpos,radius=2./100*griddims[0],color='k',alpha=1.))
	ax.add_patch(plt.Circle(minzpos,radius=2./100*griddims[0],color='k',alpha=1.))
	ax.add_patch(plt.Circle(maxzpos,radius=1./100*griddims[0],color='r',alpha=1.))
	ax.add_patch(plt.Circle(minzpos,radius=1./100*griddims[0],color='b',alpha=1.))
	if show_colorbar:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img,cax=cax)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_xlabel(r'$\mathrm{z\,(nm)\in['+str('%3.3f'%minz)+','+str('%3.3f'%maxz)+']}$',fontsize=16)
	
	ax = fig.add_subplot(gs[0,1])
	residmax = np.max(abs(residmap))
	img = ax.imshow(residmap,extent=None,origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap=mpl.cm.RdBu_r,vmin=-residmax,vmax=residmax)
	for pt in fintarget:
		ax.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='k',alpha=0.5))
	if show_colorbar:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img,cax=cax)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_title(r'$\mathrm{(z-\hat{z})\,(nm)\in['+str('%3.1f'%minresid)+','+str('%3.1f'%maxresid)+']}$',fontsize=16)

	plt.show()

#---development attempt 1 here, where I rewrite the master code block for script-preproc-dimple.py
if v1:

	from scipy import ndimage

	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

	'''
	implement PBC on domains, use extended map for the rest of this
	identify domains
	identify shadow
	program rules for selection with all geographic parameters
	plot on panels
	do the fit
	plot the fit, curvature, and domains on the same plots
	label the hmax positive and negative with text
	make a video
	add constraints for the center of the Gaussian
	maybe a bad idea to do PBC if they are all connected, so drop it. note that the PBC boundary is arbitrary
	'''

	#---NOTE: the CGMD surfaces are flush, meaning there is a duplicate row

	fr = 2
	griddims = mset.griddims
	zfilterdir = 1

	#---tesselate (via PBCs) the original surface
	#---surfpbc = numpy.tile(mset.surf[fr],(3,3))
	surfpbc = mset.surf[fr]
	#---positive/negative domain selection
	surfpoz = array([[(1 if surfpbc[i][j] > 0 else 0) 
		for j in range(1*(mset.griddims[1]-1))] for i in range(1*(mset.griddims[0]-1))]).T
	surfneg = array([[(1 if surfpbc[i][j] < 0 else 0) 
		for j in range(1*(mset.griddims[1]-1))] for i in range(1*(mset.griddims[0]-1))]).T
	label_im_poz, nb_labels = ndimage.label(surfpoz)
	label_im_neg, nb_labels = ndimage.label(surfneg)
	surf_discrete = label_im_poz - label_im_neg
	#---protein shadow selection
	prot_disc_pts = [(int(round(pt[0]/vecs[0]*(griddims[0]-1))),int(round(pt[1]/vecs[1]*(griddims[1]-1)))) 
		for pt in mset.protein[fr]]
	#---the following is equivalent
	#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	prot_locs = [[(1 if (i,j) in prot_disc_pts else 0) for i in range(mset.griddims[0]-1)] 
		for j in range(mset.griddims[1]-1)]
	#protpos = array(where(array(mset.rezipgrid(proteins[fr%len(proteins)]))!=0.0)).T
	gridpos = [[i,j] for i in range(mset.griddims[0]-1) for j in range(mset.griddims[1]-1)]
	distmat = scipy.spatial.distance.cdist(prot_disc_pts,gridpos)
	bufferlist = [gridpos[j] for j in list(set([w for v in [list(where(distmat[i]<cutoff)[0]) 
		for i in range(len(distmat))] for w in v]))]
	shadow = [[(1 if [i,j] in bufferlist else 0) for i in range(mset.griddims[0]-1)] 
		for j in range(mset.griddims[1]-1)]
	#---final selection
	target = [list(i) for i in array(where(((array(surf_discrete).T>0.)+(array(shadow).T==1)))).T]
	target = [list(i) for i in array(where((array(surf_discrete).T>0.)*(array(shadow).T==1)!=0)).T]
	if zfilterdir == 1:
		supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in target])) if i > 0]
	else:
		supplement_domains = [i for i in list(set([surf_discrete[i[0],i[1]] for i in target])) if i < 0]
	fintarget = [list([l[1],l[0]]) for l in [j for k in [list(array(where(surf_discrete==i)).T) 
		for i in supplement_domains] for j in k]]
	#---convert to xyz
	targetxyz = array([[i[0]*vecs[0]/(mset.griddims[0]-1),i[1]*vecs[1]/(mset.griddims[1]-1),
		mset.surf[fr][i[0],i[1]]] for i in fintarget])

	#---plot it
	fig = plt.figure()
	ax0 = plt.subplot2grid((1,1), (0,0))
	img = ax0.imshow(surf_discrete,extent=None,origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap=mpl.cm.jet)
	for pt in fintarget:
		ax0.add_patch(plt.Circle((pt[0],pt[1]),radius=1./2/100*griddims[0],color='r',alpha=1.))
	for pt in bufferlist:
		ax0.add_patch(plt.Circle((pt[0],pt[1]),radius=1./4/100*griddims[0],color='w',alpha=1.))
	for pt in [list(i) for i in array(where(array(prot_locs)==1)).T]:
		ax0.add_patch(plt.Circle((pt[0],pt[1]),radius=1./3/100*griddims[0],color='k',alpha=1.))
	divider = make_axes_locatable(ax0)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img,cax=cax)
	plt.show()

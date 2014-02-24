#!/usr/bin/python -i

if 1:

	from membrainrunner import *

	location = ''
	execfile('locations.py')
	execfile('plotter.py')
	
	tilt_pickles = ['pkl.tilt.membrane-v599-select.md.part0003.select.nosol.pkl',
		'pkl.tilt.membrane-v701.md.part0003.60000-160000-200.pkl',
		'pkl.tilt.membrane-v700.md.part0006.360000-460000-200.pkl',
		'pkl.tilt.membrane-v700.md.part0002.100000-200000-200.pkl']

	mset = unpickle(pickles+tilt_pickles[-1])
	
	sysname = tilt_pickles[-2][18:-4]
	
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

'''
...

'''

if 1:
	areas = mset.store[0].get(['monolayer',0,'type','area'])
	angles = mset.store[0].get(['monolayer',0,'type','angle'])
	positions = array(mset.store[1].data)[:,0]

	fig = plt.figure()

	numgridpts = 32
	vecs = mset.vecs[0]

	#---coallate data
	binned_data_angles = []
	binned_data_areas = []
	for fr in range(len(angles)):
		print fr
		combined_data = [[positions[fr][i,0],positions[fr][i,1],180.*angles[fr][0][i]] for i in range(len(angles[0][0]))]
		binned_data_angles.append(mset.rezipgrid(array(combined_data),vecs=mset.vecs[0],grid=(numgridpts,numgridpts)))
		combined_data = [[positions[fr][i,0],positions[fr][i,1],1./3.*areas[fr][0][i]] for i in range(len(angles[0][0]))]
		binned_data_areas.append(mset.rezipgrid(array(combined_data),vecs=mset.vecs[0],grid=(numgridpts,numgridpts)))

if 1:

	ax0 = plt.subplot2grid((1,2), (0,0))
	ax0.set_title('lipid tilt angle')
	maximum = max([max(i) for i in mean(binned_data_angles,axis=0)])+10
	minimum = min([min(i) for i in mean(binned_data_angles,axis=0)[1:-1][1:-1]])-10
	img = ax0.imshow(mean(binned_data_angles,axis=0), extent=None, origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap='Greys',vmax=maximum,vmin=minimum)
	for pt in np.mean(mset.protein,axis=0):
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2*1.5/64*numgridpts,color='r',
			alpha=1.0)
		ax0.add_patch(circ)
	cax = inset_axes(ax0,
         width="5%",
         height="100%",
         bbox_transform=ax0.transAxes,
         bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
         loc= 1)
	fig.colorbar(img,cax=cax)
	plt.tight_layout() 
	
	plt.savefig(pickles+'fig-'+sysname+'-lipid-tilts.png',dpi=500,bbox_inches='tight')
	
	plt.show()
	fig = plt.figure()

	ax1 = plt.subplot2grid((1,2), (0,0))
	ax1.set_title('lipid area (A2)')
	img = ax1.imshow(mean(binned_data_areas,axis=0), extent=None, origin='LowerLeft',
		interpolation='nearest',aspect='equal',cmap='Greys')
	for pt in mset.protein[0]:
		circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
			int(round(pt[1]/vecs[0]*numgridpts))),radius=1./2*1.5/64*numgridpts,color='r',
			alpha=1.0)
		ax1.add_patch(circ)

	cax = inset_axes(ax1,
         width="5%",
         height="100%",
         bbox_transform=ax1.transAxes,
         bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
         loc= 1)
	fig.colorbar(img,cax=cax)
	plt.tight_layout()
	
	plt.savefig(pickles+'fig-'+sysname+'-lipid-areas.png',dpi=500,bbox_inches='tight')
	
	plt.show()

if 1:
	'''
	outline
		open both pickles
		for each frame
			get lipids-proteins distance matrix
			get surface normals
			get 
	'''
	

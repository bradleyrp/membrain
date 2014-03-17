#!/usr/bin/python -i

from membrainrunner import *
from scipy import spatial
from scipy import linalg

location = ''
execfile('locations.py')

#---Selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
headspan = ['resname PI2P and (name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'resname PI2P and (name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)']
headangle = 'resname PI2P and (name C2 or name P or name C14)'

#---Analysis plan
analysis_descriptors = {
	'v509-40000-90000-1000':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-spanangle',
		'trajsel':'s6-kraken-md.part0018.40000-90000-1000.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Na',
		'name': 'PtdIns(4,5)P$_2$ with Na$^+$'},
	'v510-40000-75000-1000':
		{'sysname':'membrane-v510',
		'sysname_lookup':'membrane-v510-spanangle',
		'trajsel':'s8-kraken-md.part0021.40000-75000-1000.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Mg',
		'name': 'PtdIns(4,5)P$_2$ with Mg$^{2+}$'},
	'v511-40000-88000-1000':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-spanangle',
		'trajsel':'s8-kraken-md.part0021.40000-88000-1000.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Cal',
		'name': 'PtdIns(4,5)P$_2$ with Ca$^{2+}$'},
	'v533-30000-35000-100':
		{'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-spanangle',
		'trajsel':'s3-sim-kraken-md.part0011.30000-35000-100.spanangle.xtc',
		'resname':'P35P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Mg',
		'name': 'PtdIns(3,5)P$_2$ with Mg$^{2+}$'},
	'v534-50000-55000-100':
		{'sysname':'membrane-v534',
		'sysname_lookup':'membrane-v534-spanangle',
		'trajsel':'s5-trestles-md.part0017.50000-55000-100.spanangle.xtc',
		'resname':'P35P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Ca',
		'name': 'PtdIns(3,5)P$_2$ with Ca$^{2+}$'},

	}
	
		
analysis_names = ['v534-50000-55000-100']
routine = ['compute','plot','unpickle'][0:2]


#---MAIN

for aname in analysis_names:
	head_area = []
	head_angle = []
	if 'compute' in routine:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		# Needs a -pbc nojump trajectory.
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		clock = []
		result_data = MembraneData('spanangle')

		residues = mset.universe.selectAtoms('resname '+resname).resids()
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)		
			for res in residues:
				coords = mset.universe.selectAtoms('resid '+str(res)+' and resname '+resname+' and '+headspan)
				head_size = (max(spatial.distance.pdist(coords.coordinates())))
				head_area.append(head_size*head_size)
				coords = mset.universe.selectAtoms('resid '+str(res)+' and '+headangle).coordinates()
				angle = (arccos(np.dot(coords[0]-coords[1],coords[2]-coords[1])/np.linalg.norm(coords[0]-coords[1])/np.linalg.norm(coords[2]-coords[1])))
				head_angle.append(angle*(180./3.1415926))
				# For pickling.
				result_data.data.append([mset.universe.trajectory[fr].time,str(res),head_size*head_size,angle*(180./3.1415926)])
			clock.append(mset.universe.trajectory[fr].time)
		mset.store.append(result_data)
		pickle.dump(mset,open(pickles+'pkl.headspan-headangle.'+aname+'.pkl','w'))


	if 'plot' in routine:
		font = {'family' : 'sans-serif',
			'size'   : 22}
		mpl.rc('font', **font)
		mpl.rc('text', usetex=True)
		mpl.rc('text.latex', preamble='\usepackage{sfmath}')
		mpl.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}',r'\usepackage{amsmath}',
							r'\usepackage{siunitx}',r'\sisetup{detect-all}',
						        r'\usepackage{helvet}',r'\usepackage{sansmath}',
						        r'\sansmath', r'\usepackage{upgreek}']
		mpl.rcParams['xtick.major.pad'] = 8
		mpl.rcParams['ytick.major.pad'] = 8
		fig = plt.figure(figsize=(11,8.5))
		gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.05)
		ax = fig.add_subplot(gs[0])
	#	ax1.hist2d(head_angle, head_area)
	
		if 'unpickle' in routine:
			mset = unpickle(pickles+'pkl.headspan-headangle.'+aname+'.pkl')
			mdat = mset.store[0]
			head_angle.append(array(mdat.get(['headspan',0])[:,3]))
			head_area.append(array(mdat.get(['headspan',0])[:,2]))
			# This seems to work with any keyword given to mdat.get
			tmp = reduce(lambda x,y: x.extend(float(y)),head_area)
			head_area = [float(i) for i in tmp]
			tmp = reduce(lambda x,y: x.extend(float(y)),head_angle)
			head_angle = [float(i) for i in tmp]
			sysname = analysis_descriptors[aname]["sysname"]
			name = analysis_descriptors[aname]["name"]
		H, xedges, yedges = histogram2d(head_angle,head_area,bins=41,normed=True)
		midx = (xedges[1:]+xedges[:-1])/2
		midy = (yedges[1:]+yedges[:-1])/2
		extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
		cmap = mpl.cm.jet
		cmap.set_bad(cmap(0),1.)
	#	ax.imshow(array(H).T, extent=None, interpolation='nearest',aspect='equal',origin='lower',norm=None,cmap=cmap)
		ax.set_title(name)
		ax.set_xlabel('Head-tail angle (degrees)}$')
		ax.set_ylabel(r'Molecular area (\AA$^2$)')

		im = mpl.image.NonUniformImage(ax, interpolation='nearest')	
		xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
		ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
		im.set_data(xcenters, ycenters, H)
		ax.set_xlim(xedges[0], xedges[-1])
		ax.set_ylim(yedges[0], yedges[-1])
		ax.images.append(im)
		X, Y = np.meshgrid(xedges, yedges)
		plot = ax.pcolormesh(X, Y, H, figure=fig, visible=True)
		plt.colorbar(plot)
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-size-angle-correlation.png',dpi=300,bbox_inches='tight')
#		plt.show()

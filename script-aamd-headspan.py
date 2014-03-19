#!/usr/bin/python -i

from membrainrunner import *
from scipy import spatial
from scipy import linalg
import scipy.optimize as so

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
	'v509-40000-90000-100':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-spanangle',
		'trajsel':'s6-kraken-md.part0018.40000-90000-100.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Na',
		'name': 'PtdIns(4,5)P$_2$ with Na$^+$'},
	'v510-40000-75000-100':
		{'sysname':'membrane-v510',
		'sysname_lookup':'membrane-v510-spanangle',
		'trajsel':'s8-kraken-md.part0021.40000-75000-100.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Mg',
		'name': 'PtdIns(4,5)P$_2$ with Mg$^{2+}$'},
	'v511-40000-88000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-spanangle',
		'trajsel':'s8-kraken-md.part0021.40000-88000-100.spanangle.xtc',
		'resname':'PI2P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP42 or name OP43 or name OP44)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Ca',
		'name': 'PtdIns(4,5)P$_2$ with Ca$^{2+}$'},
	'v533-40000-54000-100':
		{'sysname':'membrane-v533',
		'sysname_lookup':'membrane-v533-spanangle',
		'trajsel':'s4-sim-kraken-md.part0015.40000-54000-100.spanangle.xtc',
		'resname':'P35P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Mg',
		'name': 'PtdIns(3,5)P$_2$ with Mg$^{2+}$'},
	'v534-40000-60000-100':
		{'sysname':'membrane-v534',
		'sysname_lookup':'membrane-v534-spanangle',
		'trajsel':'s4-sim-kraken-md.part0013.40000-60000-100.spanangle.xtc',
		'resname':'P35P',
		'headspan':'(name OP52 or name OP53 or name OP54 or name OP32 or name OP33 or name OP34)',
		'headangle':'(name C2 or name P or name C14)',
		'ionname':'Ca',
		'name': 'PtdIns(3,5)P$_2$ with Ca$^{2+}$'},
	}
	
		
analysis_names = ['v509-40000-90000-100','v510-40000-75000-100', 'v511-40000-88000-100', 'v533-40000-54000-100', 'v534-40000-60000-100'][0:5]
routine = ['compute','plot','unpickle','contour'][0:2]

# These functions are from: https://gist.github.com/adrn/3993992
def find_confidence_interval(x, pdf, confidence_level):
	return pdf[pdf > x].sum() - confidence_level
 
def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
	""" Create a density contour plot.
 
	Parameters
	----------
	xdata : numpy.ndarray
	ydata : numpy.ndarray
	nbins_x : int
		Number of bins along x dimension
	nbins_y : int
		Number of bins along y dimension
	ax : matplotlib.Axes (optional)
		If supplied, plot the contour to this axis. Otherwise, open a new figure
	contour_kwargs : dict
		kwargs to be passed to pyplot.contour()
	"""
 
	H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
	x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
	y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
 
	pdf = (H*(x_bin_sizes*y_bin_sizes))
 
	one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
	two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
	three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
	levels = [one_sigma, two_sigma, three_sigma]
 
	X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
	Z = pdf.T
 
	if ax == None:
		contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
	else:
		contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
 
	return contour


#---MAIN

for aname in analysis_names:
	head_area = []
	head_angle = []
	data_max = []
	if 'compute' in routine:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		# Needs a -pbc nojump trajectory.
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
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
			print "Frame "+str(fr)+" done"
		mset.store.append(result_data)
		pickle.dump(mset,open(pickles+'pkl.headspan-headangle.'+aname+'.pkl','w'))
		checktime()

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
		H, xedges, yedges = histogram2d(head_angle,head_area,bins=41,normed=True,range=((60,180),(45,100)))
		midx = (xedges[1:]+xedges[:-1])/2
		midy = (yedges[1:]+yedges[:-1])/2
		extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]
		data_max.append(H.max)
		ionname = analysis_descriptors[aname]["ionname"]
		if ionname == 'Na':
			cmap = mpl.cm.Greens
		elif ionname == 'Mg':
			cmap = mpl.cm.Reds
		elif ionname == 'Ca':
			cmap = mpl.cm.Blues
		else:
			cmap = mpl.cm.jet
#		cmap.set_bad(cmap(0),1.)
		im = ax.imshow(array(H).T, extent=(60,180,45,100), interpolation='nearest',aspect='auto',origin='lower',norm=None,cmap=cmap)
		fig.colorbar(im) 
#		im.set_clim(0,0.003)
		ax.set_title(name)
		ax.set_xlabel('Head-tail angle (degrees)}$')
		ax.set_ylabel(r'Molecular area (\AA$^2$)')
#		plt.pcolor(xedges, yedges, H) # This is totally different than imshow.
#		plt.show()
		
		plt.savefig(pickles+'fig-'+aname+'-size-angle-correlation.png',dpi=300,bbox_inches='tight')

	if 'contour' in routine:
		density_contour(head_angle, head_area, 41, 41)
		plt.show()

		

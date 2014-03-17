#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

analysis_descriptors = {

	'v514-10000-29000-100':
		{'sysname':'membrane-v514',
		'sysname_lookup':'membrane-v514-ions',
		'trajsel':'s3-sim-compbio-md.part0004.10000-29000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v514.a2-surfacer.s3-sim-compbio-md.part0004.10000-29000-100.pkl',
		'ionname':'NA'},

	'v532-20000-58000-100':
		{'sysname':'membrane-v532',
		'sysname_lookup':'membrane-v532-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-58000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',
		'ionname':'Cal'},
		
	'v531-20000-62000-100':
		{'sysname':'membrane-v531',
		'sysname_lookup':'membrane-v531-ions',
		'trajsel':'s4-sim-trestles-md.part0007.20000-62000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v531.a6-surfacer.s4-sim-trestles-md.part0007.20000-62000-100.pkl',
		'ionname':'MG'},

	'v530-30000-100000-100':
		{'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-ions',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl',
		'ionname':'NA'},
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl',
		'ionname':'Cal'}}
analysis_names = ['v531-20000-62000-100']
routine = ['compute','postproc','computexyz','plot'][3]

#---plot details
nbins = 40
allalpha = 1.
edgeprop = 2.

#---MAIN
#-------------------------------------------------------------------------------------------------------------

# This plots the ion distribution above the membranes.

if 'plot' in routine:
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	import matplotlib.gridspec as gridspec
	from matplotlib.ticker import MaxNLocator
	from scipy import optimize

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
	ax1 = fig.add_subplot(gs[0])

	hist, binedges = np.histogram(array(rewrapz[rewrapz > center]).flatten(),range=(center,center+50),bins=50)
	mid = (binedges[1:]+binedges[:-1])/2
	ax1.bar(mid-center,hist/420.,color=clrs[0], alpha=0.5, label="Na$^+$ above charged monolayer")
	
	hist, binedges = np.histogram(array(rewrapz[rewrapz < center]).flatten(),range=(center-50,center),bins=50)
	mid = (binedges[1:]+binedges[:-1])/2
	ax1.bar(center-mid,hist/420.,color=clrs[1], alpha=0.5, label="Na$^+$ above neutral monolayer")

	ax1.grid(True)
	ax1.legend(loc=1,fontsize=18)
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
	ax1.set_title('Counterion distribution',fontsize=22)
	ax1.set_ylabel('Average number',fontsize=22)
	ax1.set_xlabel(r'Distance above interpolated membrane surface (\AA)',fontsize=22)
	gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
#	plt.show()
	plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-roundz-analysis.png',dpi=300)


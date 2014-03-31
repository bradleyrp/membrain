#!/usr/bin/python

execfile('script-header.py')
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

results = [\
	'pkl.cells.membrane-v531.a6-surfacer.s4-sim-trestles-md.part0007.20000-62000-100.pkl',\
	'pkl.cells.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',\
	'pkl.cells.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl'\
	]
aname = ['v531-part0007.20000-62000-100','v532-part0007.20000-58000-100','v530-part0006.30000-100000-100']
descriptions = ["Mg$^{2+}$", "Ca$^{2+}$", "Na$^+\,$"]

for ad in range(0,3):

	fig = plt.figure(figsize=(11,8.5))
	gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.05)
	ax1 = fig.add_subplot(gs[0])

	mset = []
	mset = unpickle(pickles+results[ad])
	'''
	mset.getdata('cells').get(['monolayer',0,'type','areas'])
	available_residues = mset.resnames
	pi2p = available_residues[3]
	mset.monolayer_by_resid
	'''
	tmp = []
	tmp2 = []

	#---relative resids
	resids = [[mset.resids_reracker.index(i) for i in j] for j in mset.resids]
	#---color codes
	colorcodes = [[[i for i in range(4) if j in resids[i]][0] for j in mono] for mono in mset.monolayer_residues]
	pi2p_resids = where(np.in1d(colorcodes, 3)==True) # resids is equal to the index because it's re-racked.

	tmp = mset.getdata('cells').get(['monolayer',0,'type','areas'])
	#tmp[0][pi2p_resids[0][0]:pi2p_resids[0][-1]+1]
	tmp2 = array([tmp[i][pi2p_resids[0][0]:pi2p_resids[0][-1]+1] for i in range(len(tmp))]).flatten()
	
	hist,binedge = numpy.histogram(tmp2,bins=100,normed=True)
	mid = (binedge[1:]+binedge[:-1])/2
	ax1.plot(mid,hist,c=clrs[ad%len(clrs)],lw=4,label='PtdIns(4,5)P$_2$ with '+descriptions[ad]) # ps to ns
	ax1.fill_between(mid,hist,0,color=clrs[ad%len(clrs)],alpha=0.5)
	ax1.grid(True)
	ax1.legend(loc=0,fontsize=18)
	#plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
	#ax1.set_title('Lateral diffusion of PtdIns(4,5)P$_2$ (10\% mole fraction)',fontsize=22)
	ax1.set_xlabel('Area per molecule (\AA$^2$)',fontsize=22)
	ax1.set_ylabel('Relative frequency', fontsize=22)
	# Left align the scientific notation.
	'''
	for tick in ax1.yaxis.get_major_ticks():
	    tick.tick1line.set_markersize(0)
	    tick.tick2line.set_markersize(0)
	    tick.label1.set_horizontalalignment('left')
	yax = ax1.get_yaxis()
	yax.set_tick_params(pad=50)
	'''
	#fig.text(0, 0.5, "MSD (nm$^2$)", rotation="vertical", va="center",fontsize=22)	
	gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
	#plt.show()
	plt.savefig('/home/davids/repo-pickles/fig-'+aname[ad]+'.png',dpi=300)
	plt.close()

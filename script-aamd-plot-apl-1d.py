#!/usr/bin/python

import re

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

#results = [\
#	'pkl.cells.membrane-v531.a7-headspan.s4-sim-trestles-md.part0007.20000-62000-100.pkl',\
#	'pkl.cells.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',\
#	'pkl.cells.membrane-v530.a2-apl.u5-sim-trestles-md.part0006.30000-100000-100.pkl',
#	]
#
#results = [\
#
#	'pkl.cells.membrane-v533.a2-spanangle.s4-sim-kraken-md.part0015.40000-54000-100.pkl',\
#	'pkl.cells.membrane-v534.a1-spanangle.s4-sim-kraken-md.part0013.40000-60000-100.pkl',
#	]

results = [\
	'pkl.cells.membrane-v531.a7-headspan.s4-sim-trestles-md.part0007.20000-62000-100.pkl',\
	'pkl.cells.membrane-v532.a5-surfacer.s4-sim-trestles-md.part0007.20000-58000-100.pkl',\
	'pkl.cells.membrane-v533.a2-spanangle.s4-sim-kraken-md.part0015.40000-54000-100.pkl',\
	'pkl.cells.membrane-v534.a1-spanangle.s4-sim-kraken-md.part0013.40000-60000-100.pkl',
	]

	
#aname = ['v531-part0007.20000-62000-100','v532-part0007.20000-58000-100','v530-part0006.30000-100000-100']
aname = ['v531-part0007.20000-62000-100','v532-part0007.20000-58000-100','v533-part0015.40000-54000-100','v534-part0015.40000-54000-100']
#aname = ['v533-part0015.40000-54000-100','v534-part0015.40000-54000-100']
#descriptions = ["Mg$^{2+}$", "Ca$^{2+}$", "Na$^+$"]
descriptions = ["PtdIns(4,5)P$_2$ with Mg$^{2+}$", "PtdIns(4,5)P$_2$ with Ca$^{2+}$", "PtdIns(3,5)P$_2$ with Mg$^{2+}$", "PtdIns(3,5)P$_2$ with Ca$^{2+}$"]
#descriptions = ["Mg$^{2+}$", "Ca$^{2+}$"]

routine = ['compute','plot'][:]
keyword_re = re.compile("|".join(map(re.escape, ['v533','v534'])))

fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.05)
ax1 = fig.add_subplot(gs[0])

for ad in range(len(aname)):
	if 'compute' in routine:

		mset = []
		mset = unpickle(pickles+results[ad])
		tmp = []
		tmp2 = []

		#---relative resids
		'''
		colorcodes = [[[i for i in range(len(mset.resids)) if j in mset.resids[i]][0] for j in mono] for mono in mset.monolayer_residues]
		pi2p_resids = where(np.in1d(colorcodes[0], 4)==True) # resids is equal to the index because it's re-racked.
		pc_resids = where(np.in1d(colorcodes[1], 0)==True)
		pe_resids_lower = where(np.in1d(colorcodes[0], 2)==True)
		ps_resids = where(np.in1d(colorcodes[0], 3)==True)

		lower = mset.getdata('cells').get(['monolayer',0,'type','areas'])
		upper = mset.getdata('cells').get(['monolayer',1,'type','areas'])
		tmp2 = array([lower[i][pi2p_resids[0][0]:pi2p_resids[0][-1]+1] for i in range(len(lower))]).flatten()
		pc = array([lower[i][pc_resids[0][0]:pc_resids[0][-1]+1] for i in range(len(upper))]).flatten()
		
		pe_lower = array([lower[i][pe_resids_lower[0][0]:pe_resids_lower[0][-1]+1] for i in range(len(lower))]).flatten()
		pe_upper = []
		pe = pe_lower
		ps = array([lower[i][ps_resids[0][0]:ps_resids[0][-1]+1] for i in range(len(lower))]).flatten()
		'''
		
		
		pi2p = []
		ps = []
		pe = []
		pc = []
		chl1_bottom = []
		chl1_top = []
		mono = 0
		if not(bool(keyword_re.search(aname[ad]))):
			resn = 'PI2P'
			inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
			pi2p.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
			
		else:
			resn = 'P35P'
			inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
			pi2p.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
        	resn = 'DOPS'
		inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
        	ps.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
        	resn = 'DOPE'
		inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
        	pe.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])

		resn = 'CHL1'
		inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
		chl1_bottom.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
		mono = 1
        	resn = 'POPC'
		inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
		pc.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
		resn = 'CHL1'        	
		inds = array( [mset.monolayer_residues[mono].index(i) for i in mset.monolayer_by_resid[mono][mset.resnames.index(resn)]])
		chl1_top.append([array(i)[inds] for i in mset.getdata('cells').get(['monolayer',mono,'type','areas'])])
		
		'''
		print 'Sanity check:'
		print '# PC = '+str(shape(pc)[0]/mset.nframes)
		print '# PE (upper) = '+str(shape(pe_upper)[0]/mset.nframes)+ ' and (lower) = '+str(shape(pe_lower)[0]/mset.nframes)
		print '# PS = '+str(shape(ps)[0]/mset.nframes)
		print '# PIP2 = '+str(shape(tmp2)[0]/mset.nframes)
		'''
		print 'Sanity check:'
		print '# PC = '+str(shape(array(pc).flatten())[0]/mset.nframes)
		print '# PE = '+str(shape(array(pc).flatten())[0]/mset.nframes)
		print '# PS = '+str(shape(array(ps).flatten())[0]/mset.nframes)
		print '# PIP2 = '+str(shape(array(pi2p).flatten())[0]/mset.nframes)
		
		hist,binedge = numpy.histogram(pi2p,bins=100,normed=True,range=(20,140))
		hist_pc,binedge_pc = numpy.histogram(pc,bins=100,normed=True,range=(20,140))
		hist_pe,binedge_pe = numpy.histogram(pe,bins=100,normed=True,range=(20,140))
		hist_ps,binedge_ps = numpy.histogram(ps,bins=100,normed=True,range=(20,140))

		hist_chl1_top, binedge_chl1_top = numpy.histogram(chl1_top,bins=100,normed=True,range=(20,140))
		hist_chl1_bot, binedge_chl1_bot = numpy.histogram(chl1_bottom,bins=100,normed=True,range=(20,140))
		mid_chl1_top = (binedge_chl1_top[1:]+binedge_chl1_top[:-1])/2
		mid_chl1_bot = (binedge_chl1_bot[1:]+binedge_chl1_bot[:-1])/2

		mid = (binedge[1:]+binedge[:-1])/2
		mid_pc = (binedge_pc[1:]+binedge_pc[:-1])/2
		mid_pe = (binedge_pe[1:]+binedge_pe[:-1])/2
		mid_ps = (binedge_ps[1:]+binedge_ps[:-1])/2
	if 'plot' in routine:

		def fill_between(x, y1, y2=0, ax=None, **kwargs):
			"""Plot filled region between `y1` and `y2`.

			This function works exactly the same as matplotlib's fill_between, except
			that it also plots a proxy artist (specifically, a rectangle of 0 size)
			so that it can be added it appears on a legend.
			"""
			ax = ax if ax is not None else plt.gca()
			ax.fill_between(x, y1, y2, **kwargs)
			p = plt.Rectangle((0, 0), 0, 0, **kwargs)
			ax.add_patch(p)
			return p

		ax1.plot(mid,hist,c=clrs[ad%len(clrs)],lw=4) # ps to ns
		'''
		ax1.plot(mid_pc,hist_pc,c=clrs[3],lw=2,) # ps to ns
		ax1.plot(mid_ps,hist_ps,c=clrs[4],lw=2,) # ps to ns
		ax1.plot(mid_pe,hist_pe,c=clrs[6],lw=2,) # ps to ns
		'''
		#ax1.fill_between(mid,hist,0,color=clrs[ad%len(clrs)],alpha=0.5)
		#ax1.fill_between(mid_pc,hist_pc,0,color=clrs[3],alpha=0.3,label='asdfasdf')
		#ax1.fill_between(mid_ps,hist_ps,0,color=clrs[4],alpha=0.3)
		#ax1.fill_between(mid_pe,hist_pe,0,color=clrs[5],alpha=0.3)
		
		#ax1.set_title('Asymmetric membranes with cholesterol and '+descriptions[ad],fontsize=22)
		#fill_between(mid_pc,hist_pc,0,ax1,color=clrs[3],alpha=0.3,label='PtdCho')
		#fill_between(mid_pe,hist_pe,0,ax1,color=clrs[4],alpha=0.3,label='PtdEtn')
		#fill_between(mid_pe,hist_ps,0,ax1,color=clrs[6],alpha=0.3,label='PtdSer')
		#fill_between(mid_chl1_top,mean(array([hist_chl1_top,hist_chl1_bot]),axis=0),0,ax1,color='0.5',alpha=0.3,label='Cholesterol')
		#if bool(keyword_re.search(aname[ad])):
		#	fill_between(mid,hist,0,ax1,color=clrs[ad%len(clrs)],alpha=0.3,label='PtdIns(3,5)P$_2$')
		#else:
		#	fill_between(mid,hist,0,ax1,color=clrs[ad%len(clrs)],alpha=0.3,label='PtdIns(4,5)P$_2$')

		fill_between(mid,hist,0,ax1,color=clrs[ad%len(clrs)],alpha=0.3,label=descriptions[ad])
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
		#plt.savefig('/home/davids/repo-pickles/fig-apl-cells-'+aname[ad]+'.png',dpi=300)
		#plt.close()
#plt.show()
#plt.savefig('/home/davids/repo-pickles/fig-apl-cells-'+aname[ad]+'.png',dpi=300)
plt.savefig('/home/davids/repo-pickles/fig-apl-cells-v531-v532-v533-v534.png',dpi=300)
#plt.close()

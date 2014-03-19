#!/usr/bin/python -i

from membrainrunner import *

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = 20
framecount_force_exact = True
location = ''
execfile('locations.py')

#---imports
import matplotlib.gridspec as gridspec
#mpl.rc('text.latex', preamble='\usepackage{sfmath}')
#mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams.update({'font.size': 14})

#---colors
which_brewer_colors = [0,1,2,3,4,5,6,7]
clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]

#---analyses
analysis_descriptors = [
	('pkl.gr2d.v531.PIP2-PIP2.s4-sim-trestles-md.part0004.12000-16000-10.pkl','Mg$^{2+}$ 12-16 ns'),
	('pkl.gr2d.v531.PIP2-PIP2.s5-sim-kraken-md.part0020.54000-60000-10.pkl','Mg$^{2+}$ 54-60 ns'),
	('pkl.gr2d.v532.PIP2-PIP2.s4-sim-trestles-md.part0004.12000-16000-10.pkl','Ca$^{2+}$ 12-16 ns'),
	('pkl.gr2d.v532.PIP2-PIP2.s5-sim-kraken-md.part0018.54000-60000-10.pkl','Ca$^{2+}$ 54-60 ns'),
	('pkl.gr2d.v531.POPC-POPC.s5-sim-kraken-md.part0020.58000-62000-100.atomP.pkl','Mg$^{2+}$ 58-62 ns')
	]
#analyses = [analysis_descriptors[i] for i in range(2,4)]
analyses = [analysis_descriptors[4]]
analysis_name = 'v531-popc-popc-late'

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over analyses
allmids = []
fig = plt.figure(figsize=(8,6))
gs = mpl.gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])
for ad in analyses:
	(picklename,label) = ad
	#---load
	del mset
	mset = unpickle(pickles+picklename)
	mdat = mset.store[0]
	#---reconstitute the g(r) curves and relevant properties
	allcurves = mdat.get(['monolayer',1]) # for POPC
	binsizeabs = mdat.getnote('binsizeabs')
	'''
	allcurves = mdat.get(['monolayer',0])
	binsizeabs = mdat.getnote('binsizeabs')
	cutoff = mdat.getnote('cutoff')
	scanrange = arange(0,int(cutoff),binsizeabs)
	cutoff = max(scanrange)
	nbins=len(scanrange)-1
	'''
	cutoff_bins = min([shape(allcurves[i])[0] for i in range(len(allcurves))])
	scanrange = arange(0,int(cutoff_bins))
	nbins = len(scanrange)-1
	avgcurv = np.mean(array([i[0:nbins] for i in allcurves]),axis=0)
	hist,binedge = numpy.histogram([1 for i in range(1000)],range=(0,max(scanrange)*binsizeabs),bins=len(scanrange)-1)
	mid = (binedge[1:]+binedge[:-1])/2
	binwidth = binsizeabs
	areas = [pi*binwidth*mid[i]*2 for i in range(len(mid))]
	nlipids = mdat.getnote('points_counts')[0][0]
	vecs = np.mean(mset.vecs,axis=0)
	#---plot
	ax.plot(mid,avgcurv/areas/(nlipids/(vecs[0]*vecs[1])),'o-',c=clrs[analyses.index(ad)],label=label)
	allmids.append(mid)
ax.set_xlim((0,80))
ax.grid(True)
ax.set_xlabel(r'r (\AA)', fontsize=fsaxtext)
ax.set_ylabel(r'g(r)', fontsize=fsaxtext)
ax.set_title('PtdCho -- PtdCho pair distribution function')
plt.legend(loc='upper right')
plt.savefig(pickles+'fig-gr2d-'+analysis_name+'.png',dpi=500,bbox_inches='tight')
plt.show()


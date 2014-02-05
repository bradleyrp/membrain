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
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams.update({'font.size': 14})

#---colors
which_brewer_colors = [0,1,2,3,4,5,6,7]
clrs = [brewer2mpl.get_map('paired','qualitative',9).mpl_colors[i] for i in which_brewer_colors]

#---analyses
analysis_descriptors = [
	('pkl.gr2d.v531.PIP2-PIP2.s4-sim-trestles-md.part0004.12000-16000-10.pkl','v531 12-16 ns'),
	('pkl.gr2d.v531.PIP2-PIP2.s5-sim-kraken-md.part0020.54000-60000-10.pkl','v531 54-60 ns'),
	('pkl.gr2d.v532.PIP2-PIP2.s4-sim-trestles-md.part0004.12000-16000-10.pkl','v532 12-16 ns'),
	('pkl.gr2d.v532.PIP2-PIP2.s5-sim-kraken-md.part0018.54000-60000-10.pkl','v532 54-60 ns')
	]
analyses = [analysis_descriptors[i] for i in range(4)]
analysis_name = 'v531-v532-pi2p-early-late'

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over analyses
fig = plt.figure(figsize=(8,6))
gs = mpl.gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])
for ad in analyses:
	(picklename,label) = ad
	print label
	#---load
	del mset
	mset = unpickle(pickles+picklename)
	print len(mset.store)
	mdat = mset.store[0]
	#---reconstitute the g(r) curves and relevant properties
	allcurves = array(mdat.data)[0]
	nbins=len(arange(0,int(mdat.getnote('cutoff')),mdat.getnote('binsizeabs')))-1
	avgcurv = np.mean(array([i for i in allcurves if len(i) == nbins]),axis=0)
	print avgcurv
	cutoff = mdat.getnote('cutoff')
	nbins=len(arange(0,int(mdat.getnote('cutoff')),mdat.getnote('binsizeabs')))-1
	hist,binedge = numpy.histogram([1 for i in range(1000)],range=(0,int(cutoff)),bins=nbins)
	mid = (binedge[1:]+binedge[:-1])/2
	binwidth = mid[1]-mid[0]
	areas = [pi*binwidth*mid[i]*2 for i in range(len(mid))]
	nlipids = mdat.getnote('points_counts')[0][0]
	vecs = np.mean(mset.vecs,axis=0)
	#---plot
	#---HACKED BELOW
	ax.plot(mid,avgcurv/areas/(nlipids/(vecs[0]*vecs[1])),'o-',c=clrs[analyses.index(ad)],label=label)
ax.set_xlim((0,80))
ax.grid(True)
ax.set_xlabel(r'$\mathbf{r\:(\AA)}$',fontsize=16)
ax.set_ylabel(r'$\mathbf{g(r)}$',fontsize=16)
plt.legend(loc='lower right')
plt.savefig(pickles+'fig-'+analysis_name+'.png',dpi=500,bbox_inches='tight')
plt.show()


#!/usr/bin/python -i

from scipy import optimize
from matplotlib.ticker import MaxNLocator

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

clrs = [brewer2mpl.get_map('paired','qualitative',4).mpl_colors[i] for i in range(4)]

fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.05)
ax = fig.add_subplot(gs[0])

allrawcurves = array([mean(distsxyz[i],axis=0) for i in range(nframes)]).T
fast = []
slow = []
for c in range(len(allrawcurves)):
	curv = allrawcurves[c]
#	ax.plot(array(times)/1000.,curv,c=clrs[0] if curv[-1] > 500 else clrs[2], alpha=0.2) # Ad hoc limit.
	if curv[-1] < 500:
		ax.plot(array(times)/1000.,curv,c=clrs[0] if curv[-1] > 500 else clrs[2], alpha=0.7) # Ad hoc limit.
	fast.append(curv) if curv[-1] > 500 else slow.append(curv)
#ax.plot(array(times)/1000.,mean(allrawcurves,axis=0),c='k', alpha=1, label='Mean of all data')
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y: (y - fitfunc(p, x))
powerlaw = lambda x, amp, index: amp * (x**index)
pinit = [100.0, -100.0]

fast_out = optimize.leastsq(errfunc, pinit, args=(log10(times[1:]), log10(mean(fast[:][1:],axis=0)[1:])), full_output=1)
slow_out = optimize.leastsq(errfunc, pinit, args=(log10(times[1:]), log10(mean(slow[:][1:],axis=0)[1:])), full_output=1)
#ax.plot(array(times)/1000., mean(fast[:][1:],axis=0), c='k', lw=2, label='Mean fast')
#ax.plot(array(times)/1000., mean(slow[:][1:],axis=0), c='k', lw=2, label='Mean slow')

pfinal = fast_out[0]
covar = fast_out[1]
index = pfinal[1]
amp = 10.0**pfinal[0]
indexErr = sqrt( covar[0][0] )
ampErr = sqrt( covar[1][1] ) * amp
print index, indexErr, amp, ampErr
#ax.plot(array(times)/1000., powerlaw(times, amp, index), c=clrs[1], lw=4, label='Power law fit, amplitude = %.2f, exponent = %.2f' %(amp,index))
pfinal = slow_out[0]
covar = slow_out[1]
index = pfinal[1]
amp = 10.0**pfinal[0]
print index, amp
ax.plot(array(times)/1000., powerlaw(times, amp, index), c=clrs[3], lw=4, label='Power law fit, amplitude = %.2f, exponent = %.2f' %(amp,index))
ax.plot(array(times)/1000., powerlaw(times, amp, float64(1.0)), c='k', lw=4, label='Power law fit, amplitude = %.2f (forced), exponent = %.2f' %(amp,index))
#ax.set_ylim((10**-1,10**6))
ax.grid(True)
ax.legend(loc=2,fontsize=18)
plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
#ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlabel('Time (ns)',fontsize=22)
# Left align the scientific notation.
for tick in ax.yaxis.get_major_ticks():
    tick.tick1line.set_markersize(0)
    tick.tick2line.set_markersize(0)
    tick.label1.set_horizontalalignment('left')
yax = ax.get_yaxis()
yax.set_tick_params(pad=50)
fig.text(0, 0.5, "MSD (nm$^2$)", rotation="vertical", va="center",fontsize=22)	
gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
#plt.show()
#plt.savefig('/home/davids/Science/p1-bilayers/Analysis/Ion-diffusion/v511-MSD-slow-forcedfit.png',dpi=300)
plt.clf()
plt.close()

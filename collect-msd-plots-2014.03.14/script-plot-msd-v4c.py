#!/usr/bin/python2.7 -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 6).mpl_colors
#clrs = brewer2mpl.get_map('paired', 'qualitative', 6).mpl_colors
import scipy.stats as stats
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

base = "../membrane-v5xx/"
# Rearranging so that calcium is blue, magnesium red, sodium green.
#systems = ["membrane-v509/a3-lipidwise-msd/", "membrane-v510/a2-msd-for-real/", "membrane-v511/a1-msd/"]
systems = [ "membrane-v510/a2-msd-for-real/", "membrane-v511/a1-msd/", "membrane-v509/a3-lipidwise-msd/",]
descriptions = ["Mg$^{2+}$", "Ca$^{2+}$", "Na$^+\,$"]
pi2p_analyses = ["msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg"]
dops_analyses = ["msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg"]
dopcpopc_analyses = ["msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg"]

fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(1,1,wspace=0.0,hspace=0.05)
ax1 = fig.add_subplot(gs[0])

# This variant does not assume a linear fit and instead reports sub/super diffusion.

for ad in range(0,3):
	# PIP2
	lines = []
	fp = open(base+systems[ad]+pi2p_analyses[ad],'r')
	D = ''
	for line in fp:
		if line[0] != '#' and line[0] != '@':
			lines.append([float(j) for j in line.split()])
		if len(line) > 5: 
			if line[2] == 'D' and line[3] == '[':
				D_gromacs = line[19:25]
	fp.close()
	lines = array(lines)
	tdatx = array(lines)[1:,0]
	tdaty = array(lines)[1:,1]
	first_ten = int(len(lines)*.1)+2
	last_ten =  int(len(lines)*.9)+1
	logx = log10(tdatx[first_ten:last_ten])
	logy = log10(tdaty[first_ten:last_ten])
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	powerlaw = lambda x, amp, index: amp * (x**index)
	pinit = [100.0, -100.0]
	out = optimize.leastsq(errfunc, pinit, args=(logx, logy), full_output=1)
	pfinal = out[0]
 	covar = out[1]
	print pfinal
	print covar
	index = pfinal[1]
	amp = 10.0**pfinal[0]
 	indexErr = sqrt( covar[0][0] )
	ampErr = sqrt( covar[1][1] ) * amp
	print amp
	print index
	ax1.plot(tdatx/1000., powerlaw(tdatx, amp, index), c='k', lw=2)
	
	
#	[q,r], covariance = polyfit(tdatx[first_ten:last_ten], tdaty[first_ten:last_ten], 1, cov=True)
#	print 'Gromacs D = ' + str(float(D_gromacs)*10**3) + ' (micron**2/second)' + ' Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
#	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
	# Plot
	ax1.plot(tdatx/1000.,tdaty,'-',c=clrs[ad%len(clrs)],lw=2,alpha=1.0,zorder=10, label=descriptions[ad]+', Amplitude = %.2f, Exponent = %.2f' % (1000*amp, index)) # ps to ns
	ax1.grid(True)
	ax1.legend(loc=2,fontsize=18)
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
	ax1.set_title('Lateral diffusion of PtdIns(4,5)P$_2$ (10\% mole fraction)',fontsize=22)
	ax1.set_yscale('log')
	ax1.set_xscale('log')
#	ax1.set_xticklabels([])
	ax1.set_xlabel('Time (ns)',fontsize=22)
	# Left align the scientific notation.
	for tick in ax1.yaxis.get_major_ticks():
	    tick.tick1line.set_markersize(0)
	    tick.tick2line.set_markersize(0)
	    tick.label1.set_horizontalalignment('left')
	yax = ax1.get_yaxis()
	yax.set_tick_params(pad=50)
	
fig.text(0, 0.5, "MSD (nm$^2$)", rotation="vertical", va="center",fontsize=22)	
gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
#plt.show()
plt.savefig('/home/davids/repo-pickles/fig-v509v510v511-powerlaw.png',dpi=300)
plt.close()

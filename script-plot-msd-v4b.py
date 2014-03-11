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
#descriptions = ["10\% PIP2, 72% DOPC, Na$^+$", "Symmetric, 10\% PIP2, Mg$^{2+}$", "Symmetric, 10\% PIP2, Ca$^{2+}$"]
descriptions = ["Mg$^{2+}$", "Ca$^{2+}$", "Na$^+\,$"]
pi2p_analyses = ["msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg"]
dops_analyses = ["msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg"]
dopcpopc_analyses = ["msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg"]

fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(3,1,wspace=0.0,hspace=0.05)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

# This variant plots all the PI2P diffusion curves on separate panels for each ion.
# With a black line for slope (linear fit!).

for sys in range(0,3):
	ax = plt.subplot(gs[sys])
	all_this_x = []
	all_this_y = []
	for ad in range(0,80):
		# PIP2
		lines = []
		name="msd-lateral-pi2p-%04d.xvg" %ad
		fp = open(base+systems[sys]+name,'r')
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
		all_this_x.append([lines[first_ten:last_ten,0][i] for i in range(len(lines[first_ten:last_ten]))])
		all_this_y.append([lines[first_ten:last_ten,1][i] for i in range(len(lines[first_ten:last_ten]))])
#		print 'Gromacs D = ' + str(float(D_gromacs)*10**3) + ' (micron**2/second)' + ' Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
			# Plot
		ax.plot(tdatx/1000.,tdaty,'-',c=clrs[sys%len(clrs)],lw=2,alpha=0.2,zorder=1, label=descriptions[sys] if ad == 0 else "__nolegend__")
#		ax.plot(tdatx[first_ten:last_ten]/1000.,fitted_y, c='k', lw=2)
		ax.grid(True)
		ax.legend(loc=2,fontsize=18)
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
#		ax1.set_title('Lateral diffusion of PtdIns(4,5)P$_2$ (10\% mole fraction)',fontsize=22)
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_ylim(10**-3, 10**2)
		# Left align the scientific notation.
		for tick in ax.yaxis.get_major_ticks():
		    tick.tick1line.set_markersize(0)
		    tick.tick2line.set_markersize(0)
		    tick.label1.set_horizontalalignment('left')
		yax = ax.get_yaxis()
		yax.set_tick_params(pad=50)
		if sys != 2:
			ax.set_xticklabels([])
		else:
			ax.set_xlabel('Time (ns)',fontsize=22)
	flat_x = [x for sublist in all_this_x for x in sublist]
	flat_y = [y for sublist in all_this_y for y in sublist]
	[q,r], covariance = polyfit(flat_x[1:], flat_y[1:], 1, cov=True)
	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
	ax.plot(tdatx[first_ten:last_ten]/1000.,fitted_y, c='k', lw=2, alpha=1, label="D = %.2f $\upmu \mathsf{m}^2$/s" % float(q*10**6/4.))
	# Make lines in the legend solid (not alpha = whatever).
	leg = plt.legend(loc=2, prop={'size':18})
	for l in leg.get_lines():
	    l.set_alpha(1)
	

fig.text(0, 0.5, "MSD (nm$^2$)", rotation="vertical", va="center",fontsize=22)	
gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
plt.show()
#plt.savefig('/home/davids/Science/p1-bilayers/Analysis/Lipid-diffusion/PI2P-by-ion.png',dpi=300)
#plt.close()

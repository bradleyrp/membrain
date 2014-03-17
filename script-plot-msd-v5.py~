#!/usr/bin/python2.7 -i

import matplotlib as mpl
mpl.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 6).mpl_colors
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
systems = [ "membrane-v510/a2-msd-for-real/", "membrane-v511/a1-msd/", "membrane-v509/a3-lipidwise-msd/",]
descriptions = ["Mg$^{2+}$", "Ca$^{2+}$", "Na$^+\,$"]
pi2p_analyses = ["msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg"]
dops_analyses = ["msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg"]
dopcpopc_analyses = ["msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg"]

fig = plt.figure(figsize=(11,8.5))
#gs = gridspec.GridSpec(3,1,wspace=0.0,hspace=0.2)
#ax1 = fig.add_subplot(gs[0])
#ax2 = fig.add_subplot(gs[1])
#ax3 = fig.add_subplot(gs[2])
gs = gridspec.GridSpec(1,1)
ax = plt.subplot(gs[0])

# This variant plots the PI2P diffusion curves on top of each other.

last_delta = []
for sys in range(0,3):
#	ax = plt.subplot(gs[sys])
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
		last_delta.append(tdaty[-1])
		first_ten = int(len(lines)*.1)+2
		last_ten =  int(len(lines)*.9)+1
		# Plot the whole thing, but only fit 10-90% of the total data.
		[q,r], covariance = polyfit(tdatx[first_ten:last_ten], tdaty[first_ten:last_ten], 1, cov=True)
		print 'Gromacs D = ' + str(float(D_gromacs)*10**3) + ' (micron**2/second)' + ' Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
		fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
		# Plot
		#ax1.plot(tdatx/1000.,tdaty,'-',c=clrs[ad%len(clrs)],lw=2,alpha=1.0,zorder=10, label=descriptions[ad]+", D = %.2f $\upmu \mathsf{m}^2$/s" % float(q*10**6/4.))
		ax.plot(tdatx/1000.,tdaty,'-',c=clrs[sys%len(clrs)],lw=2,alpha=0.2,zorder=10)
		ax.grid(True)
		ax.legend(loc=2,fontsize=18)
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
#		ax1.set_title('Lateral diffusion of PtdIns(4,5)P$_2$ (10\% mole fraction)',fontsize=22)
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_ylim(10**-3, 10**2)
#		if sys != 2:
#			ax.set_xticklabels([])


fig.text(0, 0.5, "MSD (nm$^2$)", rotation="vertical", va="center",fontsize=22)	
ax.set_xlabel('Time (ns)',fontsize=22)
#gs.tight_layout(fig, rect=[0.03, 0.0, 1, 1]) # Leave space for the common y-label.
plt.show()
#plt.savefig('/tmp/tmp6.png',dpi=300)
#plt.close()
#plt.figure()
#plt.hist(log10(last_delta), color = 'blue') 
#plt.savefig('/tmp/tmp7.png',dpi=300)
#plt.close()

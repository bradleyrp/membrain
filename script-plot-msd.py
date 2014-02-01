#!/usr/bin/python -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rc('text.latex', preamble='\usepackage{siunitx}')
mpl.rc('text', usetex=True)

# The following command is broken for me.
#mpl.rcParams.update({'font.style':'sans-serif'})
mpl.rcParams.update({'font.size': 16})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

import scipy.stats as stats


#basedir = '/home/davids/xfer/v532/' # Change this later.
basedir = '/home/rpb/worker/lipidwise-msd/'
analysis_descriptors = [str(basedir)+'msd-'+str("%04d"%i)+'.xvg' for i in range(0,40)] 
analyses = [analysis_descriptors[i] for i in range(len(analysis_descriptors))]

fig = plt.figure(figsize=(10,4))
gs = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
axes = []
ax = fig.add_subplot(gs[0,0])
axes.append(ax)
bzs = []
all_y = []
all_x = []
for ad in analyses:
	#---load
	lines = []
	fp = open(ad,'r')
	for line in fp:
		if line[0] != '#' and line[0] != '@':
			lines.append([float(i) for i in line.split()])
	fp.close()
	#---plot
	lines = array(lines)
	tdatx = array(lines)[1:,0]
	tdaty = array(lines)[1:,1]
	# Create a list of all the y values (MSD).
	first_ten = int(len(lines)*.1)+2
	last_ten =  int(len(lines)*.9)+1
	all_x.append([lines[first_ten:last_ten,0][i] for i in range(len(lines[first_ten:last_ten]))])
	all_y.append([lines[first_ten:last_ten,1][i] for i in range(len(lines[first_ten:last_ten]))])
	tdatall = log10(lines[1:])
	tdat = log10(lines[int(len(tdatall)*.1)+2:int(len(tdatall)*.9)+1])
	[bz,az] = polyfit(tdat[:,0],tdat[:,1],1)
	print str(ad)+' D = '+str(bz)
	bzs.append(bz)
	fitted = [bz*i+az for i in tdatall[:,0]]
	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,alpha=0.5)
	if 0:
		ax.plot(tdat[:,0],tdat[:,1],'-',c='r',lw=2)
		ax.plot(tdatall[:,0],fitted,'-',c='k',lw=1,alpha=0.5)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
ax.set_xscale('log')
ax.set_yscale('log')
#--- Recreate the best fit for all the curves, create a 95% confidence interval and shade.
all_x_flat = [log10(x) for sublist in all_x for x in sublist]
all_y_flat = [log10(y) for sublist in all_y for y in sublist]

[q,r], covariance = polyfit(all_x_flat, all_y_flat, 1, cov=True)
print 'Overall average D = ' + str(q) + ' (some crazy units)'

fitted_y = [q*i+r for i in log10(lines[first_ten:last_ten,0])]
fitted_x = lines[first_ten:last_ten,0]
# Confidence interval:
t = stats.t.ppf(0.975, len(fitted_y) - 2) # Students' t distribution 97.5 percentile, n-2 d.o.f.
# The problem is we need to find the residual for all the y values to the fitted line, but we can't do a point-by-point subtraction because fitted is smaller than all_y_flat. Also, it won't work to average all y values for a single x point and do that subtraction, because that is losing degrees of freedom and thus, power.
#residuals = [all_y_flat[i]] - fitted[i] for i in range(len(fitted))]
#s_err = sqrt(sum(residuals**2)/(len(all_y_flat) - 2))  # Standard deviation of the error (residuals)
#ci = t * s_err * sqrt(1/len(all_y_flat) + (x2 - mean(all_x_flat))**2/sum((all_x_flat-mean(all_x_flat))**2))
unlogged_y = [10**fitted_y[i] for i in range(len(fitted_y))]
ax.plot(lines[first_ten:last_ten,0],unlogged_y,'-',c='k',lw=2,label='D = %1.4f'%q)
plt.legend()
#ax.fill_between(unlogged_x,unlogged_y+ci,unlogged_y-ci, c='k',alpha=0.1)


#---histogram
ax = fig.add_subplot(gs[1])
axes.append(ax)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_xscale('linear')
ax.set_yscale('linear')

nbins = 10
histo,binedge = histogram(bzs,bins=nbins,range=(0,max(bzs)*1.1), density=True)
# Setting density=True does not totally normalize the histogram, unless the bin width = 1.
mid = (binedge[1:]+binedge[:-1])/2.
# Set align='center' to center the bins on the midpoint.
ax.bar(mid,histo/sum(histo),color=clrs[0],alpha=0.7,width=max(bzs)*1.1/nbins, align='center')
# Label every other bin center.
ax.xaxis.set_ticks(mid[1:-1:2])
# Units are totally weird here. nm^2/ps = 10^{6}*micron^2/second
# These D should be around 1-20 micron^2/second
# What we have are a million times bigger.
ax.set_xlabel(r'D (\SI{1e6}{\micro\metre^2\per\second})')

plt.show()



#!/usr/bin/python2.7 -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
#plt.rc('font', **{'sans-serif': 'Arial'})

import scipy.stats as stats

base = "../membrane-v5xx/"
systems = ["membrane-v509/a3-lipidwise-msd/", "membrane-v510/a2-msd-for-real/", "membrane-v511/a1-msd/", "membrane-v530/a1-msd/", "membrane-v531/a4-full-on-msd/", "membrane-v532/a3-full-on-msd/"]
pi2p_analysis = ["msd-lateral-pi2p-all.xvg"]
dops_analysis = ["msd-lateral-dops-all.xvg"]
dopcpopc_analysis = ["msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-popc-all.xvg", "msd-lateral-popc-all.xvg", "msd-lateral-popc-all.xvg"]

#analysis_strings = ["PI2P (80 molecules)", "DOPC (576 molecules)"]
#analyses = [systems[i] for i in range(len(systems))]
fig, ax = plt.subplots(nrows=1)

for ad in systems:
	#---load
	lines = []
	fp = open(base+ad+dopcpopc_analysis[systems.index(ad)],'r')
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
	# Plot the whole thing, but only fit 10-90% of the total data.
	[q,r], covariance = polyfit(tdatx[first_ten:last_ten], tdaty[first_ten:last_ten], 1, cov=True)
	print 'Overall average D = ' + str(q*10**6/4.) + ' (micron**2/second)'
	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
#	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,alpha=1.0, label=(analysis_strings[analyses.index(ad)%len(analysis_strings)])+ ', D = %s $\mu$m$^2\!$/s' % round( (q*10**6/4.),2), zorder=10)
	ax.plot(tdatx,tdaty,'-',c=clrs[systems.index(ad)%len(clrs)],lw=2,alpha=1.0,zorder=10)
#	ax.plot(tdatx[first_ten:last_ten], fitted_y, '-', c='k', lw=3, alpha=1, zorder=1) 
	
ax.set_xlabel('Time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
plt.legend(loc=2, prop={'size':18})
plt.show()



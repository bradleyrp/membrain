#!/usr/bin/python2.7 -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
#clrs = brewer2mpl.get_map('Set1', 'qualitative', 6).mpl_colors
clrs = brewer2mpl.get_map('Dark2', 'qualitative', 6).mpl_colors

fig_width_pt = 1024
fig_height_pt = 1024
inches_per_pt = 1.0/72.27               # Convert pt to inches
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_height_pt*inches_per_pt       # height in inches
fig_size = [fig_width,fig_height]
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Sans']})
#params = {'backend': 'ps',
#          'text.latex.preamble': [r"\usepackage{sansmath}",
#          			  r"\usepackage{helvet}",
#          			  r"\usepackage{siunitx}",
#          			  r"\sansmath",
#         			  r"\renewcommand{\familydefault}{\sfdefault}",
#         			  r"\sisetup{math-rm=\mathsf,text-rm=\sffamily}"],
#          'axes.labelsize': 24,
#          'text.fontsize': 24,
#          'legend.fontsize': 18,
#          'xtick.labelsize': 20,
#          'ytick.labelsize': 20,
#          'text.usetex': True,
#          'figure.figsize': fig_size,
#          'axes.unicode_minus': True}
#plt.rcParams.update(params) 


import scipy.stats as stats

base = "../membrane-v5xx/"
# Main 6 now.
systems = ["membrane-v509/a3-lipidwise-msd/", "membrane-v510/a2-msd-for-real/", "membrane-v511/a1-msd/", "membrane-v530/a1-msd/", "membrane-v531/a4-full-on-msd/", "membrane-v532/a3-full-on-msd/"]
descriptions = ["Symmetric, 10% PIP2, Na$^+$", "Symmetric, 10% PIP2, Mg$^{2+}$", "Symmetric, 10% PIP2, Ca$^{2+}$", "Asymmetric, 5% PIP2, Na$^+$", "Asymmetric, 5% PIP2, Mg$^{2+}$", "Asymmetric, 5% PIP2, Ca$^{2+}$"]
pi2p_analyses = ["msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", "msd-lateral-pi2p-all.xvg", ]
dops_analyses = ["msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", "msd-lateral-dops-all.xvg", ]
dopcpopc_analyses = ["msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-dopc-all.xvg", "msd-lateral-popc-all.xvg", "msd-lateral-popc-all.xvg", "msd-lateral-popc-all.xvg"]

fig, ax = plt.subplots(nrows=1)

for ad in range(0,6):
	#---load
	lines = []
	fp = open(base+systems[ad]+pi2p_analyses[ad],'r')
	D = ''
	for line in fp:
		if line[0] != '#' and line[0] != '@':
			lines.append([float(j) for j in line.split()])
		if len(line) > 5: 
			if line[2] == 'D' and line[3] == '[':
				D = line[19:25]
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
	print 'Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
#	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,alpha=1.0, label=(analysis_strings[analyses.index(ad)%len(analysis_strings)])+ ', D = %s $\mu$m$^2\!$/s' % round( (q*10**6/4.),2), zorder=10)
	ax.plot(tdatx,tdaty,'-',c=clrs[ad%len(clrs)],lw=2,alpha=1.0,zorder=10, label=descriptions[ad]+", D = "+str(float(D)*1000.)+"$\mathsf{\mu m} ^2$/second")
#	ax.plot(tdatx[first_ten:last_ten], fitted_y, '-', c='k', lw=3, alpha=1, zorder=1) 
	
	
ax.set_xlabel('Time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
plt.legend(loc=2, prop={'size':18})
plt.show()

fig, ax = plt.subplots(nrows=1)
for ad in range(0,6):
	#---load
	lines = []
	fp = open(base+systems[ad]+dops_analyses[ad],'r')
	D = ''
	for line in fp:
		if line[0] != '#' and line[0] != '@':
			lines.append([float(j) for j in line.split()])
		if len(line) > 5: 
			if line[2] == 'D' and line[3] == '[':
				D = line[19:25]
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
	print 'Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
#	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,alpha=1.0, label=(analysis_strings[analyses.index(ad)%len(analysis_strings)])+ ', D = %s $\mu$m$^2\!$/s' % round( (q*10**6/4.),2), zorder=10)
	ax.plot(tdatx,tdaty,'-',c=clrs[ad%len(clrs)],lw=2,alpha=1.0,zorder=10, label=descriptions[ad]+", D = "+str(float(D)*1000.)+" $\mu m^2$/second")
#	ax.plot(tdatx[first_ten:last_ten], fitted_y, '-', c='k', lw=3, alpha=1, zorder=1) 
	
	
ax.set_xlabel('Time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
plt.legend(loc=2, prop={'size':18})
plt.show()

fig, ax = plt.subplots(nrows=1)
for ad in range(0,6):
	#---load
	lines = []
	fp = open(base+systems[ad]+dopcpopc_analyses[ad],'r')
	D = ''
	for line in fp:
		if line[0] != '#' and line[0] != '@':
			lines.append([float(j) for j in line.split()])
		if len(line) > 5: 
			if line[2] == 'D' and line[3] == '[':
				D = line[19:25]
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
	print 'Fitted D = ' + str(q*10**6/4.) + ' (micron**2/second)'
	fitted_y = [q*i+r for i in tdatx[first_ten:last_ten]]
#	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,alpha=1.0, label=(analysis_strings[analyses.index(ad)%len(analysis_strings)])+ ', D = %s $\mu$m$^2\!$/s' % round( (q*10**6/4.),2), zorder=10)
	ax.plot(tdatx,tdaty,'-',c=clrs[ad%len(clrs)],lw=2,alpha=1.0,zorder=10, label=descriptions[ad]+", D = "+str(float(D)*1000.)+" $\mu$m$^2$/second")
#	ax.plot(tdatx[first_ten:last_ten], fitted_y, '-', c='k', lw=3, alpha=1, zorder=1) 
	
	
ax.set_xlabel('Time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
plt.legend(loc=2, prop={'size':18})
plt.show()
#!/usr/bin/python -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import brewer2mpl
clrs = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams.update({'font.style':'sans-serif'})
mpl.rcParams.update({'font.size': 16})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

basedir = '/home/rpb/queue/lipidwise-msd/' # Change this later.
analysis_descriptors = [str(basedir)+'msd-'+str("%04d"%i)+'.xvg' for i in range(0,40)] 
analyses = [analysis_descriptors[i] for i in range(len(analysis_descriptors))]
analyses = analysis_descriptors[:10]

fig = plt.figure(figsize=(10,4))
gs = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0)
axes = []
ax = fig.add_subplot(gs[0])
axes.append(ax)
bzs = []
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
	tdatall = log10(lines[1:])
	tdat = log10(lines[int(len(tdatall)*.1)+2:int(len(tdatall)*.9)+1])
	[bz,az] = polyfit(tdat[:,0],tdat[:,1],1)
	print bz
	bzs.append(10**-5*bz)
	fitted = [bz*i+az for i in tdatall[:,0]]
	ax.plot(tdatx,tdaty,'-',c=clrs[analyses.index(ad)%len(clrs)],lw=2,label=ad,alpha=0.8)
	if 0:
		ax.plot(tdat[:,0],tdat[:,1],'-',c='r',lw=2)
		ax.plot(tdatall[:,0],fitted,'-',c='k',lw=1,alpha=0.5)
ax.set_xlabel('time (ps)')
ax.set_ylabel('MSD (nm$^{2}$)')
ax.set_xscale('log')
ax.set_yscale('log')

#---histogram
ax = fig.add_subplot(gs[1])
axes.append(ax)
nbins = 10
ax.hist(bzs,nbins,normed=1,facecolor=clrs[0],alpha=0.8)
'''
histo,binedge = histogram(bzs,bins=nbins,range=(0,max(bzs)*1.1))
mid = (binedge[1:]+binedge[:-1])/2.
ax.bar(mid,histo,color=clrs[0],alpha=0.7,width=max(bzs)*1.1/nbins)
ax.set_ylim((0,max(histo)*1.1))
ax.yaxis.tick_right()
ax.set_ylabel('frequency',fontsize=18,rotation=270)
ax.yaxis.set_label_position("right")
#ax.fill_between(mid,[0 for i in histo],histo,facecolor='b',alpha=0.25,where=histo>[0 for i in histo])
ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
'''
plt.show()



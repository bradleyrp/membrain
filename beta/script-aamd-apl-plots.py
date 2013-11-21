#!/usr/bin/python -i

if 1:
	from membrainrunner import *
	execfile('plotter.py')

#---Load
#-------------------------------------------------------------------------------------------------------------

if 0:
	v509a = unpickle('pkl.apl-cells.membrane-v509.md.part0010.pkl')
	v509b = unpickle('pkl.apl-cells.membrane-v509.md.part0035.pkl')
	v514 = unpickle('pkl.apl-cells.membrane-v514.md.part0008.pkl')
	
if 1:
	v509 = unpickle('pkl.apl-cells.membrane-v509.md.part0035.pkl')
	v510 = unpickle('pkl.apl-cells.membrane-v510.md.part0033.pkl')
	v511 = unpickle('pkl.apl-cells.membrane-v511.md.part0033.pkl')

#---Load
#-------------------------------------------------------------------------------------------------------------

def plot_apl_lipids(mset,sysid=None,ltypes=None,colorset=None,plotby='frame'):
	'''Individual plots of lipid areas.'''
	if ltypes == None:
		ltypes = mset.resnames
	for typ in ltypes:
		print ltypes
		typno = mset.resnames.index(typ)
		areas = mset.getdata('cells').get(['type','areas'])
		areas3 = [[[areas[fr][mono][mset.monolayer_residues[mono].index(i)] for i in mset.monolayer_by_resid[mono][typno]] for mono in [0,1]] for fr in range(len(areas))]
		if plotby == 'monolayer':
			plotdat = mean(areas3,axis=(0,2)) 
		elif plotby == 'lipid':
			plotdat = mean(areas3,axis=(0,1))
		elif plotby == 'frame':
			plotdat = mean(areas3,axis=(1,2)) 
		elif plotby == 'all':
			plotdat = [t for u in [i for j in areas3 for i in j] for t in u]
		print shape(plotdat)
		sysid = sysid + ' ' if sysid != None else ''
		thisplot = hist(plotdat,alpha=0.5,color=colorset[typno],
			label=sysid+typ+'\n'+r"$A_{"+str(typ)+"}= "+str('%2.2f'%mean(plotdat))+
			"\pm"+str('%1.2f'%std(plotdat))+"$")
	plt.legend(loc=1)
	plt.show()
	
def plot_apl_lipids_stack(msets,sysids=None,ltypes=None,colorset=None,plotby='frame',xlims=None,normed=True,coloroffset=0):
	'''Multiple stacked plots of lipid areas.'''
	numplots = len(msets)
	for m in range(len(msets)):
		mset = msets[m]
		sysid = sysids[m]
		ax = plt.subplot2grid((numplots,1), (m,0))
		if ltypes == None:
			ltypes = mset.resnames
		for typ in ltypes:
			typno = mset.resnames.index(typ)
			areas = mset.getdata('cells').get(['type','areas'])
			areas3 = [[[areas[fr][mono][mset.monolayer_residues[mono].index(i)] for i in mset.monolayer_by_resid[mono][typno]] for mono in [0,1]] for fr in range(len(areas))]
			if plotby == 'monolayer':
				plotdat = mean(areas3,axis=(0,2)) 
			elif plotby == 'lipid':
				plotdat = mean(areas3,axis=(0,1))
			elif plotby == 'frame':
				plotdat = mean(areas3,axis=(1,2)) 
			elif plotby == 'all':
				plotdat = [t for u in [i for j in areas3 for i in j] for t in u]
			sysid = sysid + ' ' if sysid != None else ''
			# "-1" below for ion colors
			thisplot = hist(plotdat,alpha=0.5,color=colorset[(typno+(m-typno-1)*coloroffset)%len(mset.resnames)],normed=normed,
				label=sysid+typ+'\n'+r"$A_{"+str(typ)+"}= "+str('%2.2f'%mean(plotdat))+
				"\pm"+str('%1.2f'%std(plotdat))+"$")
		if xlims != None:
			ax.set_xlim(xlims)
		if m == 0:
			print 'm is zero sometime'
			plt.ylabel('Normed frequency',labelpad = 10)
		if m == numplots-1:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_position(('outward', 20))
			ax.spines['left'].set_position(('outward', 30))
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('bottom')
			plt.xlabel('Area per lipid ($\AA^{2}$)',labelpad = 10,fontsize=16)
		else:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['left'].set_position(('outward', 30))
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('bottom')
		plt.legend(loc=2)
	plt.show()

#---Main
#-------------------------------------------------------------------------------------------------------------


#---Parameters for plotting
#msets = [v509a,v509b]
#sysnames = ['early','late']
msets = [v509,v510,v511]
sysnames = ['Na','Mg','Ca']
colorset = brewer2mpl.get_map('Set1','qualitative',3).mpl_colors

#---All individual plots
if 0:
	for m in range(len(msets)):
		mset = msets[m]
		sysid = sysnames[m]
		plot_apl_lipids(mset,sysid=sysid,ltypes=None,colorset=colorset)

#---All plots, stacked
if 1:
	plot_apl_lipids_stack(
		msets,sysids=sysnames,
		ltypes=['PI2P'],
		colorset=colorset,
		xlims=[53,68],
		normed=True,
		plotby='frame',
		coloroffset=1)


########################################

#---Individual plots, per-lipid
if 0:
	mset = v509a
	sysid = sysnames[0]
	ltypes = mset.resnames
	for typ in ltypes:
		typno = mset.resnames.index(typ)
		areas = mset.getdata('cells').get(['type','areas'])
		areas3 = [[[areas[fr][m][mset.monolayer_residues[m].index(i)] for i in mset.monolayer_by_resid[0][typno]] for m in [0]] for fr in range(len(areas))]
		plotby = 'lipid'
		if plotby == 'monolayer':
			plotdat = mean(areas3,axis=(0,2)) 
		elif plotby == 'lipid':
			plotdat = mean(areas3,axis=(0,1)) 
		elif plotby == 'frame':
			plotdat = mean(areas3,axis=(1,2)) 
		elif plotby == 'all':
			plotdat = [t for u in [i for j in areas3 for i in j] for t in u]
		thisplot = hist(plotdat,normed=True,alpha=0.5,color=clrs[typno],
			label=(sysid+' ' if sysid != None else None)+typ+'\n'+r"$A_{"+str(typ)+"}= "+str('%2.2f'%mean(plotdat))+
			"\pm"+str('%1.2f'%std(plotdat))+"$"
			)
	plt.legend()
	plt.show()

if 0:
	#plot_aplxxxx(self,label='triangles',ltypes=None,specs=None,filename=None,savefig=False,monos=[0,1]):
	'''Area per lipid calculations and plot.'''
	if type(label) != str:
		label = self.avail()[label]
	lipid_areas = []
	if ltypes == None:
		ltypes = range(len(self.resnames))
	for ltype in ltypes:
		areas = []
		for mono in monos:
			print 'Calculating areas: lipid = '+self.resnames[ltype]+' monolayer '+str(mono)+'.'
			for fr in range(len(self.getdata('triangles').label)):
				points = array(self.getdata(label).get([['type','points'],['frame',fr],['monolayer',mono]]))
				tri1 = array(self.getdata(label).get([['type','lines'],['frame',fr],['monolayer',0]]))
				validres = [self.monolayer_residues[mono].index(i) 
					for i in self.monolayer_by_resid[mono][ltype]]
				areas.append([sum([1./6*norm(cross(t[0][0:2]-t[1][0:2],t[2][0:2]-t[1][0:2])) 
					for t in [[points[j] for j in tri1[i]] for i in list(where(any(tri1==p,axis=1))[0])]]) 
					for p in validres])
		lipid_areas.append(areas)
	areas = array(lipid_areas)
	fig, ax = plt.subplots()
	histplots = []
	for d in range(len(lipid_areas)):
		c = clrs[d%len(clrs)]
		thisplot = hist(mean(lipid_areas[d],axis=0),normed=True,label=self.resnames[d],alpha=0.5,color=c)
		histplots.append(thisplot)
		ax.text(0.7,0.6-d*0.05, r"$A_{"+str(self.resnames[d])+"}= "+str('%2.2f'%mean(lipid_areas[d]))+
			"\pm"+str('%1.2f'%std(lipid_areas[d]))+"$",transform=ax.transAxes,fontsize=12)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.xlabel('Area per lipid ($\AA^{2}$)',labelpad = 10,fontsize=16)
	plt.ylabel('Normed frequency',labelpad = 10,fontsize=16)
	plt.title('Area per lipid',fontsize=16)
	fig.tight_layout()
	plt.legend()
	if savefig: 
		if not filename:
			plt.savefig(name+'.png', dpi=100)
		else:
			plt.savefig(filename+'.png', dpi=100)
	plt.show()
	
	




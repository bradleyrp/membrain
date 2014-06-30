#!/usr/bin/python

showplots = False
norm_z_concentration = True
resname_group = None

#---concentration calculations

def concentration_profiler(disctraj,binedges,mset,normed=False,buffer_size=20,raw=False):
	'''Basic function for returning a concentration profile in the correct units.'''
	#---histogram
	counts,edges = histogram(
		disctraj.flatten(),
		range=(0,len(binedges[0])-1),
		bins=len(binedges[0])-1)
	counts = array(counts)/float(mset.nframes)
	meanedges = array(mean(array(binedges),axis=0))
	mids = 1./2*(meanedges[1:]+meanedges[:-1])
	vecs = array([mset.vec(i) for i in range(mset.nframes)])
	#---nions = conc mol/L * 1000L / m3 * m3 / (10**10)**3 A3 * Av particles/M * vol A3
	concfac = (6.0221413*10**23)/((10**10)**3)*10**3
	boxvol = product(mean(vecs,axis=0))
	conc = sum(counts)/boxvol/concfac
	status('status: ['+str(ion_name)+'] = '+'{:.3f}'.format(conc*1000).ljust(8)+'mM')
	conc_mM = counts/concfac/product(mean(vecs,axis=0)[:2])/binwidth*1000
	#---buffer in Angstroms at maximum distance from bilayer to compute bulk
	meanbins = mean(binedges,axis=0)
	topbin = where(meanbins>(meanbins.max()-buffer_size))[0][0]
	botbin = where(meanbins<(-meanbins.max()+buffer_size))[0][-1]
	bulkconc = mean([mean(counts[1:botbin]),mean(counts[topbin:-1])])
	bulkconc_mM = bulkconc/concfac/product(mean(vecs,axis=0)[:2])/binwidth*1000
	status('status: ['+str(ion_name)+'],bulk = '+'{:.3f}'.format(bulkconc_mM).ljust(8)+'mM')
	#---override to output raw counts
	if raw: return mids,counts
	else: return mids,(conc_mM/bulkconc_mM if normed else conc_mM)

#---concentration modular plotting functions

def plotter_concentration(mids,concs,ax=None,ion_name=None,color=None):
	'''Plot the concentration in mM units.'''
	ax.plot(mids,conc_mM,c=color,lw=2,label=proper_ion_labels[ion_name])
	ax.set_ylabel(r'$\mathrm{concentration}\:\mathrm{(mM)}$',fontsize=fsaxlabel)
	ax.set_xlabel(r'$z\:(\mathrm{\AA})$',fontsize=fsaxlabel)
	ax.grid(True)
def plotter_concentration_relative(mids,concs,mset,ax=None,ion_name=None,color=None):
	'''Plot the concentration in mM units.'''
	ax.plot(mids,concs,c=color,lw=2,label=proper_ion_labels[ion_name])
	ax.set_ylabel('relative concentration',fontsize=fsaxlabel)
	ax.set_xlabel(r'$z\:(\mathrm{\AA})$',fontsize=fsaxlabel)
	ax.grid(True)
def plotter_concentration_counts(disctraj,binedges,mset,ax=None,ion_name=None,color=None):
	'''Plot the concentration relative to the bulk.'''
	#---perform the computation here
	mids,counts = concentration_profiler(disctraj,binedges,mset,raw=True)
	centerbin = where(abs(mids)<10**-6)[0][0]
	ax.plot(mids[centerbin:],ion_charges[inum]*cumsum(counts[centerbin:]),
		c=color,lw=2,label=proper_ion_labels[ion_name])
	ax.plot(mids[:centerbin],ion_charges[inum]*cumsum(counts[:centerbin][::-1])[::-1],c=color,lw=2)
	ax.set_ylabel('cumulative counts',fontsize=fsaxlabel)
	ax.set_xlabel(r'$z\:(\mathrm{\AA})$',fontsize=fsaxlabel)
	ax.grid(True)



#---choose a system
aname = 'v532-40000-90000-50'
anum = analysis_names.index(aname)
mset = msets[anum]
	
#---prepare figure	
fig = plt.figure(figsize=(8,8))
axes = []
gs = gridspec.GridSpec(3,1)
for i in range(3):
	ax = fig.add_subplot(gs[i])
	axes.append(ax)

#---prepare figure settings
ion_colors = ['r','b']
ion_charges = [1,1]
ion_resnames = [ion_name,ion_name_alt]

#---loop over ions and make plots of concentration
for inum in range(len(ionlist)):
	#---copmpute discrete (binwise) ion trajectories
	relpos = bilayer_shift(mset=mset,vecs=array([mset.vec(i) for i in range(mset.nframes)]),
		ionspos=master_ionspos[anum][inum],midz=master_midz[anum])
	disctraj = discretizer_z(binedges,array(relpos)[valid_inds])	
	#---plot concentrations three different ways
	mids,conc_mM = concentration_profiler(disctraj,binedges,mset,normed=False)
	plotter_concentration(mids,conc_mM,ax=axes[0],ion_name=ion_resnames[inum],color=ion_colors[inum])
	mids,conc_rel = concentration_profiler(disctraj,binedges,mset,normed=True)
	plotter_concentration_relative(mids,conc_rel,mset,
		ax=axes[1],ion_name=ion_resnames[inum],color=ion_colors[inum])
	plotter_concentration_counts(disctraj,binedges,mset,
		ax=axes[2],ion_name=ion_resnames[inum],color=ion_colors[inum])
#---plot settings
for ax in axes: ax.legend(loc='center left',bbox_to_anchor=(1, 0.5),ncol=1,fancybox=True,shadow=True)
gs.tight_layout(fig,h_pad=0.6,w_pad=0.6)
#---save
plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-ion_equilibrium_z-'+aname+\
	'.png',dpi=300,bbox_inches='tight')
if showplots: plt.show()
plt.close(fig)

	


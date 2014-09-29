#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

#---inset axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

analysis_names = ['v531-40000-90000-50-bridge','v532-40000-90000-50-bridge']
target_atomsel_via_regex = 'O'
splitter = 4
routine = ['plot_residence_time_distn']

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def neighborfinder(cd,num):
	'''
	Note that by convention the "subject" is the atom referred to by the second dimension.
	'''
	natoms2 = shape(cd)[1]
	#---only do the argsort once
	aa = argsort(cd,axis=0)
	#---get the nearest lipid atoms in each nearest distance rank
	ndists = zeros((num,natoms2))
	nnexts = zeros((num,natoms2))
	for nnum in range(num):
		nnexts[nnum] = aa[nnum]
		ndists[nnum] = cd[nnexts[nnum].astype(int)][arange(natoms2),arange(natoms2)]
	return ndists,nnexts.astype(int)

def nearby(nlist,splitter):
	'''
	If distinct then only return next-ranked lipid residue numbers if they belong to a unique neighbor.
	'''
	nears = [[] for n in range(len(nlist))]
	for nnum in range(len(nlist)): nears[nnum] = where(nlist[nnum]<splitter)[0]
	return array(nears)

def ontorus(pts1,pts2,vecs):
	'''Compute distances under periodic boundary conditions with the topology of a torus.'''
	cdx = scipy.spatial.distance.cdist(pts1[...,0:1],pts2[...,0:1])
	cdy = scipy.spatial.distance.cdist(pts1[...,1:2],pts2[...,1:2])
	cdz = scipy.spatial.distance.cdist(pts1[...,2:],pts2[...,2:])
	cdx = cdx-1*(cdx>vecs[0]/2.)*vecs[0]
	cdy = cdy-1*(cdy>vecs[1]/2.)*vecs[1]
	cdz = cdz-1*(cdz>vecs[2]/2.)*vecs[2]
	cd = sqrt(cdx**2+cdy**2+cdz**2)
	return cd

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'mset' not in globals():
	all_bridged,all_contact_atoms,all_contact_lipids = [],[],[]
	all_mids,all_counts,all_durations = [],[],[]
	#---loop over systems for comparison	
	for aname in analysis_names:
		ion_name = analysis_descriptors[aname]['ion_name']
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load the lipid atoms
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = MembraneSet()
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		mset.identify_monolayers(director)
		#---load the ions
		if residues == 'infer': residues = list(set(mset.universe.residues.resnames()))
		mset.identify_residues(residues)
		mset_ions = MembraneSet()
		grofile, trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			keytrajsel='ions_trajsel',keysysname='ions_sysname')
		mset_ions.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		group_ion = mset_ions.universe.selectAtoms('name ' + str(ion_name))
		mset_ions.selections.append(group_ion)
		targlist = list(set(mset.universe.selectAtoms('all').names()))
		#---simple way to select targets based on regex for a letter in the atom name
		#---? may need replaced with more specific method
		target_names = ' or '.join(['name '+i for i in targlist if re.search(target_atomsel_via_regex,i)])
		group_lipid = mset.universe.selectAtoms(target_names)
		mset.selections.append(group_lipid)
		#---? why does switching the mset to the correct option actually change the result ??!??
		#---? check correct order
		atomnames = list(set(mset.selections[0].names()))
		#---loop over frames
		st = time.time()
		neardists = zeros((mset.nframes,len(group_ion)))
		nlatoms = len(atomnames)
		bridged,contact_atoms,contact_lipids,dists = [],[],[],[]
		for fr in range(mset.nframes):
			status('fr = '+str(fr),start=st,i=fr,looplen=mset.nframes)
			pts1 = array(mset.get_points(fr,selection_index=0))
			pts2 = array(mset_ions.get_points(fr,selection_index=0))
			#---pairwise distances between lipid atoms and ions
			cd = ontorus(pts1,pts2,mset.vec(fr))
			#---distances and atom number for the top two nearest neighbors
			ndists,nnexts = neighborfinder(cd,2)
			neardists[fr] = ndists[0]
			#---ion numbers for those within the cutoff
			nears = nearby(ndists,splitter)
			#---lipid residue numbers for nearby lipids
			nres = [nnexts[i][nears[1]]/nlatoms for i in range(2)]
			#---ions which bridge two distinct lipids
			bridged.append(nears[1][where(nres[0]!=nres[1])[0]])
			#---atoms which these ions contact
			contact_atoms.append((nnexts[:2,bridged[-1]]%nlatoms).T)
			#---contacting lipids residues
			contact_lipids.append((nnexts[:2,bridged[-1]]/nlatoms).T)
		all_bridged.append(bridged)
		all_contact_atoms.append(contact_atoms)
		all_contact_lipids.append(contact_lipids)
		#---compute residence times
		nearest = neardists.T
		durs = []
		for resnr in range(shape(nearest)[0]):
			exits = (1*(nearest[resnr][1:]<splitter)==1)+1*((nearest[resnr][:-1]<splitter)==0)
			exit_times = array(where(exits==2)[0])
			durs.append(exit_times[1:]-exit_times[:-1])
		all_durations.append(durs)
		#---histogram residence times
		lims = [0,max([max(i) for i in durs if len(i)>0])+1]
		counts = []
		for fr in range(len(durs)):
			status('fr = '+str(fr),looplen=mset.nframes)
			c,edges = histogram(durs[fr],range=lims,bins=lims[1])
			counts.append(c)
		mids = (edges[1:]+edges[:-1])/2.
		all_mids.append(mids)
		all_counts.append(counts)

#---combined residence time plot
if 'plot_residence_time_distn' in routine and 0:
	cutoff = where(mids<500)[0][-1]
	ax = plt.subplot(111)
	ax.set_title(r'$\:'+'{:.1f}'.format(splitter)+'\:\mathrm{\AA}$'+' to Phosphorus residence times',
		fontsize=fsaxlabel)
	for anum in range(len(analysis_names)):
		counts = all_counts[anum]
		mids = all_mids[anum]
		moment = mean(mean(counts,axis=0)[:cutoff]*mids[:cutoff])/mean(mean(counts,axis=0)[:cutoff])
		aname = analysis_names[anum]
		ion_name = analysis_descriptors[aname]['ion_name']
		color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
		label = analysis_descriptors[aname]['composition_name']+', '+analysis_descriptors[aname]['ion_label']+\
			'\n'+'{:.2f}'.format(moment)+' ps'
		ax.plot(mids[:20],mean(counts,axis=0)[:20],color=color,lw=3,
			label=label)
	ax.legend(loc='upper right',fontsize=fsaxlabel)
	ax.set_ylabel('frequency',fontsize=fsaxlabel)
	ax.set_xlabel('time (ps)',fontsize=fsaxlabel)
	ax.grid(True)
	plt.savefig(pickles+'fig-residence-'+'-'.join(analysis_names)+'.png',dpi=300)
	plt.show()

ax = plt.subplot(111)
ax.set_title(r'$\:'+'{:.1f}'.format(splitter)+'\:\mathrm{\AA}$'+' to Phosphorus residence times',
	fontsize=fsaxlabel)
for anum in range(len(analysis_names)):
	aname = analysis_names[anum]
	ion_name = analysis_descriptors[aname]['ion_name']
	nbridged = array([len(i) for i in all_bridged[anum]])
	c,edges = histogram(nbridged,range=(0,nbridged.max()),bins=nbridged.max())
	mids = (edges[1:]+edges[:-1])/2.
	color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
	label = analysis_descriptors[aname]['composition_name']+', '+analysis_descriptors[aname]['ion_label']
	ax.plot(mids,c,label=label,color=color,lw=3)
	ax.legend(loc='upper right',fontsize=fsaxlabel)
	ax.set_ylabel('frequency',fontsize=fsaxlabel)
	ax.set_xlabel('number of bridged ions',fontsize=fsaxlabel)
	ax.grid(True)
plt.savefig(pickles+'fig-bridge_number-'+'-'.join(analysis_names)+'.png',dpi=300)
plt.show()


if 0:
	#---plot bridged-or-not trajectories
	cbridged = concatenate(all_bridged[0])
	hbridge = histogram(concatenate(bridged),range=(0,cbridged.max()),bins=cbridged.max())
	any_bridged_ions = list(set(cbridged))
	bridge_trajectories = [[1 if inum in j else 0 for j in all_bridged[0]] for inum in any_bridged_ions]
	for inum in range(len(any_bridged_ions))[:2]:
		plt.plot(range(mset.nframes),bridge_trajectories[inum],alpha=0.2)
	plt.show()




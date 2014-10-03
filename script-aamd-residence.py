#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

#---inset axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---explicit arrays
if 0: numpy.set_printoptions(threshold='nan')

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

analysis_names = [
	'v531-40000-90000-50-bridge',
	'v532-40000-90000-50-bridge',
	][:]
target_atomsel_via_regex = 'O'
target_lipid_resname = 'PI2P'
splitter = 3.5
routine = [
	'plot_residence_time_distn',
	'plot_number_bridged_distn',
	'plot_example_bridged_trajectories',
	][:-1]

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
	all_bridged,all_contact_atoms,all_contact_lipids,all_lipid_lipid = [],[],[],[]
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
		target_names = ' or '.join(['name '+i for i in targlist if re.search(target_atomsel_via_regex,i)])
		group_lipid = mset.universe.selectAtoms('resname '+target_lipid_resname+' and ('+target_names+')')
		mset.selections.append(group_lipid)
		atomnames = list(set(mset.selections[0].names()))
		#---loop over frames
		st = time.time()
		neardists = zeros((mset.nframes,len(group_ion)))
		nlatoms,niatoms = len(atomnames),len(group_ion)
		bridged,contact_atoms,contact_lipids,dists,lipid_lipid = [],[],[],[],[]
		for fr in range(mset.nframes):
			status('fr = '+str(fr),start=st,i=fr,looplen=mset.nframes)
			pts1 = array(mset.get_points(fr,selection_index=0))
			pts2 = array(mset_ions.get_points(fr,selection_index=0))
			#---pairwise distances between lipid atoms and ions
			cd = ontorus(pts1,pts2,mset.vec(fr))
			#---distances and atom number for the top two nearest neighbors
			#---the first dimension of ndists,nnexts is the rank, the second is either the distance or atom
			ndists,nnexts = neighborfinder(cd,nlatoms+1)
			#---save the first nearest distances for later
			neardists[fr] = array(ndists[0])
			#---new method
			ndists,nnexts = ndists.T,nnexts.T
			#---see if subsequent residues are different than the first nearest residue
			subseq = array([nnexts[:,0]/nlatoms!=nnexts[:,i]/nlatoms for i in range(1,nlatoms+1)]).T
			#---since we have included pairwise distances up to nlatoms+1, each row has a true value
			#---note that this may be a bottleneck
			truenext = array([where(i==True)[0][0] for i in subseq])+1
			bridged.append(where((1*(ndists[:,0]<splitter)+\
				1*(ndists[arange(niatoms),truenext]<splitter))==2)[0])
			contact_lipids.append(
				array([(nnexts[arange(niatoms),0]/nlatoms)[bridged[-1]],
				(nnexts[arange(niatoms),truenext]/nlatoms)[bridged[-1]]]).T)
			contact_atoms.append(
				array([(nnexts[arange(niatoms),0]%nlatoms)[bridged[-1]],
				(nnexts[arange(niatoms),truenext]%nlatoms)[bridged[-1]]]).T)
			#---search for lipid-lipid contacts
			cdll = ontorus(pts1,pts1,mset.vec(fr))
			ndistsll,nnextsll = neighborfinder(cdll,nlatoms+1)
			ndistsll,nnextsll = ndistsll.T,nnextsll.T
			subseq = array([nnextsll[:,0]/nlatoms!=nnextsll[:,i]/nlatoms for i in range(1,nlatoms+1)]).T
			truenextll = array([where(i==True)[0][0] for i in subseq])+1
			closes = where(ndistsll[arange(len(ndistsll)),truenextll]<splitter*2)[0]
			lls = array([nnextsll[arange(len(ndistsll)),truenextll][closes]/nlatoms,closes/nlatoms]).T
			lipid_lipid.append(list(set([tuple(i) for i in lls if i[0]<i[1]])))
				
		#---consolidate
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
		
#---PLOT
#-------------------------------------------------------------------------------------------------------------

#---combined residence time plot
if 'plot_residence_time_distn' in routine:
	cutoff = where(mids<500)[0][-1]
	ax = plt.subplot(111)
	ax.set_title(r'$\:'+'{:.1f}'.format(splitter)+'\:\mathrm{\AA}$'+' radius residence times',
		fontsize=fsaxlabel)
	for anum in range(len(analysis_names)):
		counts = all_counts[anum]
		mids = all_mids[anum]
		moment = mean(mean(counts,axis=0)[:cutoff]*mids[:cutoff])/mean(mean(counts,axis=0)[:cutoff])
		aname = analysis_names[anum]
		ion_name = analysis_descriptors[aname]['ion_name']
		color = color_dictionary_aamd(ionname=ion_name,comparison='ions')
		label = analysis_descriptors[aname]['composition_name']+', '+\
			analysis_descriptors[aname]['ion_label']+\
			'\n'+'{:.2f}'.format(moment)+' ps'
		ax.plot(mids[:20],mean(counts,axis=0)[:20],color=color,lw=3,
			label=label)
	ax.legend(loc='upper right',fontsize=fsaxlabel)
	ax.set_ylabel('frequency',fontsize=fsaxlabel)
	ax.set_xlabel('time (ps)',fontsize=fsaxlabel)
	ax.grid(True)
	plt.savefig(pickles+'fig-residence-'+'-'.join(analysis_names)+'.png',dpi=300)
	plt.show()

#---plot distribution of bridged ion counts
if 'plot_number_bridged_distn' in routine:
	ax = plt.subplot(111)
	ax.set_title(r'$\:'+'{:.1f}'.format(splitter)+'\:\mathrm{\AA}$'+' radius residence times',
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

#---plot bridged-or-not trajectories
if 'plot_example_bridged_trajectories' in routine:
	cbridged = concatenate(all_bridged[0])
	hbridge = histogram(concatenate(bridged),range=(0,cbridged.max()),bins=cbridged.max())
	any_bridged_ions = list(set(cbridged))
	bridge_trajectories = [[1 if inum in j else 0 for j in all_bridged[0]] for inum in any_bridged_ions]
	for inum in range(len(any_bridged_ions))[:2]:
		plt.plot(range(mset.nframes),bridge_trajectories[inum],alpha=0.2)
	plt.show()


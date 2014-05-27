#!/usr/bin/python

logfile,interact,debugmode = [None,False,None]
from membrainrunner import *
execfile('locations.py')
execfile('header-aamd.py')

#---DEFAULTS
#-------------------------------------------------------------------------------------------------------------

if 'batch_override' not in globals():

	#---parameters
	rounder = 4

	#---settings
	compare_phosphate_position = False

	#---standard selection
	analysis_names = [
		'v530-40000-90000-50',
		'v531-40000-90000-50',
		'v532-40000-90000-50',
		'v533-40000-90000-50',
		'v534-40000-90000-50',
		'v509-40000-90000-50',
		'v510-40000-90000-50',
		'v511-40000-90000-50',
		'v514-22000-32000-10',
		'v515-20000-30000-10',
		][:5]
	
	#---alternate tests
	if compare_phosphate_position:
		analysis_names = [
			'v531-40000-90000-50',
			'v532-40000-90000-50',
			'v533-40000-90000-50',
			'v534-40000-90000-50',
			]
	
	#---routine
	routine = [
		'calc',
		'load',
		'plot',
		][:1]
	
	#---settings
	showplots = False
	
	#---selector
	selector = ' and '.join(list(set([key_atom_selector[i] 
		for i in key_atom_selector.keys() if i != 'CHL1'])))


#---CALCULATE
#-------------------------------------------------------------------------------------------------------------

if routine == str: routine = [routine]

if 'calc' in routine:
	#---loop over analysis questions
	for aname in analysis_names:
		#---details
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		#---loop over trajectory files
		for traj in trajfile:
			isfile = unpickle(pickles+'pkl.structures.'+'space'+str(rounder)+'A.'+\
				specname_pickle(sysname,traj)+'.pkl')
			if isfile == True:
				status('status: pkl already exists ')
			else:
				mset = MembraneSet()
				#---load the trajectory
				status('status: accessing '+specname_pickle(sysname,traj))
				mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
				checktime()
				#---average structure calculation
				mset.identify_monolayers(director)
				#---infer the timeslice from the XTC name if not specified
				if 'timeslice' in analysis_descriptors[aname].keys():
					status('warning: requested timeslice '+str(timeslice))
				else:
					if len(trajsel.split('.')[-2].split('-')) == 1:
						tslicepos = -3
						subset_name = trajsel.split('.')[-2]
					else: tslicepos = -2
					timeslice = [int(i) for i in trajsel.split('.')[tslicepos].split('-')]
				if protein_select == None:
					mset.midplaner(selector,skip=skip,rounder=float(rounder),
						framecount=framecount,timeslice=timeslice,
						thick=True)
				else:
					mset.midplaner(selector,skip=skip,rounder=float(rounder),framecount=framecount,
						protein_selection=protein_select,timeslice=timeslice)
				mset.calculate_undulations()
				#---save the data
				pickledump(mset,'pkl.structures.'+'space'+str(rounder)+'A.'+\
					specname_pickle(sysname,traj)+'.pkl',directory=pickles)
				if erase_when_finished:
					del mset
				checktime()

#---LOAD
#-------------------------------------------------------------------------------------------------------------

if 'calc' not in routine and 'msets' not in globals():
	msets = []
	#---loop over analysis questions
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset = unpickle(pickles+'pkl.structures.'+'space'+str(rounder)+'A.'+\
			specname_pickle(sysname,trajfile[0])+'.pkl')
		msets.append(mset)

#---PLOT
#-------------------------------------------------------------------------------------------------------------

if 'plot' in routine:
	fig = plt.figure()
	gs = gridspec.GridSpec(1,1)
	ax = plt.subplot(gs[0])
	for m in range(len(msets)):
		aname = analysis_names[m]
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		mset = msets[m]
		#mset.calculate_undulations(peri=True,
		#	removeavg=(True if batch != 'phosphate_position' else False))
		mset.calculate_undulations(peri=False,
			removeavg=True)
		resname = ptdins_resname
		color = color_dictionary_aamd(ionname=ion_name,lipid_resname=resname,
			comparison='ions')
		if compare_phosphate_position:
			color = color_dictionary_aamd(ionname=ion_name,lipid_resname=resname,
				comparison='ions_phospate_position')
		if batch == 'protonation':
			color = color_dictionary_aamd(lipid_resname=resname,
				comparison='protonation')
		plotter_undulate(mset,qmagfilter=[0.3,1],ax=ax,inset2d=False,
			colorspec=color,showkappa=False,
			label=ptdins_label+', '+ion_label,
			peri=peristalsis,
			intrinsic=[1,2],
			lims_override=((0.2,13),(2*10**-7,0.02)))
	ax.set_title(composition_name+' bilayer',fontsize=fsaxtitle)
	ax.legend(loc='upper right')
	plt.savefig(pickles+'/PREPARE-FIGURES/'+'fig-'+('bending-' if not peristalsis else 'peristalsis-')+\
		'-'.join(analysis_names)+'.png',\
		dpi=300,bbox_inches='tight')
	if showplots: plt.show()
	plt.close(fig)



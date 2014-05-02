#!/usr/bin/python

#---EXPERIMENT PARAMETERS, DIMPLE
#-------------------------------------------------------------------------------------------------------------

#---choose midplane resolution
gridspacing = 1.0
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'

'''
Experiment/analysis parameters list
	1. decay_z0 : 
		zero = dimple decays to zero (average height) at infinity
		min = dimple decays to the average bilayer minimum height at infinity
	2. cutoff : cutoff about protein points which defines a neighborhood (in Angstroms)
	3. z_filter : 
		inclusive = include all points
		up = include only above average (z>0) points
		down = include only below average (z<0) points
	4. fit_dynamics_type : 
		dynamic = fit a dimple to each frame
		mean = fit a dimple to the average structure
	5. framecounts : if None use all frames, if slice or integer, use that to get a subset of frames
	6. geography : 
		proteins = use the proteins themselves as geographic sources
		lattice = sample random points and use a fixed single point as the center of the neighborhood
'''

#---experiment parameters
params_expt = {
	'decay_z0' : 'zero',
	'cutoff' : 100,
	'fit_dynamics_type' : 'dynamic',
	'z_filter' : 'inclusive',
	'framecounts' : 500, 
	'geography' : '',
	}
	
#---sweeping parameters
params_sweeps = {
	'cutoff':[100],
	'geography':[
		'proteins',
		['control','v614-120000-220000-200'],
		['control','v612-75000-175000-200'],
		'max','min',
		'max_dynamic','min_dynamic',
		['lattice_grid',6],
		][:],
	}	
	
#---PLOT PARAMETERS, DIMPLE
#-------------------------------------------------------------------------------------------------------------

#---general parameters
show_plots = True

#---general settings
params_plot_settings = {
	'hifilt' : 0.06,
	'smallfilt' : 0.001,
	'hist_step' : 0.005,
	'extent_range' : 32,
	'filter_type' : 'std',
	}

#---needs notes on filter types

#---master list of parameters
params_plot_master,params_plot_master_names = [],[]

#---plot parameter set "lattice6"
params_plot_master_names.append('lattice6')
params_plot_master.append([
	[(aname,
	dict({
		'decay_z0' : 'zero',
		'cutoff' : 100,
		'fit_dynamics_type' : 'dynamic',
		'z_filter' : 'inclusive',
		'framecounts' : 500, 
		'geography' : geog,
		},))
		for geog in [['lattice_grid',6]]]
	for aname in analysis_names])

#---plot parameter set "proteins_controls"
params_plot_master_names.append('proteins_controls')
params_plot_master.append(
	[[(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			},))
		for geog in ['proteins']] 
			for aname in analysis_names 
			if (analysis_descriptors[aname])['nprots'] > 0]+\
	[[(aname,
	dict({
		'decay_z0' : 'zero',
		'cutoff' : 100,
		'fit_dynamics_type' : 'dynamic',
		'z_filter' : 'inclusive',
		'framecounts' : 500, 
		'geography' : geog,
		},))
	for geog in [['control','v614-120000-220000-200'],['control','v612-75000-175000-200']]
		for aname in analysis_names
		if (analysis_descriptors[aname])['nprots'] == 0]])

#---plot parameter set "extrema"
params_plot_master_names.append('extrema')
params_plot_master.append(
	[[(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			},))
		for geog in ['max_dynamic','min_dynamic','max','min']] 
			for aname in analysis_names])


#---plot parameter set "radar_proteins"			
params_plot_master_names.append('radar_proteins')
params_plot_master.append(
	[[(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			'neighborhood_name' : \
				('oligomer' if (analysis_descriptors[aname])['nprots'] > 1 else 'monomer'),
			},))
		for geog in ['proteins']] 
			for aname in analysis_names 
			if (analysis_descriptors[aname])['nprots'] > 0]+\
	[[(aname,
	dict({
		'decay_z0' : 'zero',
		'cutoff' : 100,
		'fit_dynamics_type' : 'dynamic',
		'z_filter' : 'inclusive',
		'framecounts' : 500, 
		'geography' : geog,
		},))]
	for geog in [['control','v614-120000-220000-200'],['control','v612-75000-175000-200']]
		for aname in analysis_names
		if (analysis_descriptors[aname])['nprots'] == 0])

#---plot parameter set "radar_extrema"			
params_plot_master_names.append('radar_extrema')
params_plot_master.append(
	[[(aname,
		dict({
			'decay_z0' : 'zero',
			'cutoff' : 100,
			'fit_dynamics_type' : 'dynamic',
			'z_filter' : 'inclusive',
			'framecounts' : 500, 
			'geography' : geog,
			},))
		for geog in ['max_dynamic','min_dynamic','max','min']] 
			for aname in analysis_names])

#!/usr/bin/python

#---EXPERIMENT PARAMETERS, TOPOGRAPHY
#-------------------------------------------------------------------------------------------------------------

#---choose midplane resolution
gridspacing = 1.0
spacetag = 'space'+str(int(round(gridspacing*10,0)))+'A.'

#---experiment parameters
params_expt = {
	'cutoff' : 10,
	}
	
#---sweeping parameters
params_sweeps = {
	'cutoff':[10],
	'geography':[
		'proteins',
		['control','v614-120000-220000-200'],
		['control','v612-75000-175000-200'],
		'max','min',
		'max_dynamic','min_dynamic',
		['lattice_grid',6],
		][3:6+1],
	}	
	
#---PLOT PARAMETERS, TOPOGRAPHY
#-------------------------------------------------------------------------------------------------------------

#---general parameters
show_plots = True

#---general settings
params_plot_settings = {
	'hifilt' : 0.06,
	'smallfilt' : 0.001,
	'hist_step' : 0.005,
	'extent_range' : 32,
	}

#---master list of parameters
params_plot_master,params_plot_master_names = [],[]

#---plot parameter set "extrema"
params_plot_master_names.append('extrema')
params_plot_master.append(
	[[(aname,
		dict({
			'cutoff' : 10,
			'geography' : geog,
			},))
		for geog in ['max_dynamic','min_dynamic','max','min']] 
			for aname in analysis_names])
			
#---snapshots for combined topography plot
viewtype = 'birdseye'
imagelist = {
	'v616-210000-310000-200':'snap-v616-'+viewtype+'.png',
	'v614-120000-220000-200':'snap-v614-'+viewtype+'.png',
	'v612-75000-175000-200':'snap-v612-'+viewtype+'.png',
	'v550-400000-500000-160':'snap-v550-'+viewtype+'.png'
	}
viewtype = 'side'
imagelist2 = {
	'v616-210000-310000-200':'snap-v616-'+viewtype+'.png',
	'v614-120000-220000-200':'snap-v614-'+viewtype+'.png',
	'v612-75000-175000-200':'snap-v612-'+viewtype+'.png',
	'v550-400000-500000-160':'snap-v550-'+viewtype+'.png'
	}


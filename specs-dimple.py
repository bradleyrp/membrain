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

'''
DIMPLE-FITTING PROCEDURE

There are currently three possible filters implemented in the 'plot_dimple' plotting section.
Here is a workflow description of the filtering.

how to read the file names:
	the "v550" numbers indicate which simulation i.e. 616 is 8xENTH, 550 is control, etc
	the cut200 indicates a buffer or cutoff of 200 Angstroms around the neighborhood
	the "mod1" indicates the filter type
		std = report maximum mean curvature in the neighborhood+buffer
		mod1 = H_max comes from the center of the dimple if it is within 1/2 of the box vector magnitude
		mod2 = H_max is the center of the dimple, but only if it's in the neighborhood+buffer

1. Select the set of interpolated bilayer points which will be fit.
	A. Choose the location of the points (the neighborhood).
		(1) The protein oligomer or monomer, imported from another simulation on the control
		(2) A single point selected from an evenly-spaced lattice.
		(3) The average or instantaneous peak or valley of the midplane
		Note: these selections are made with PBCs, so if the peak is on the edge, you still get a contiguous neighborhood.
	B. Choose a buffer around the location
	C. The fitted points will be the points inside of region defined by the neighborhood + buffer
		Note: this means that the lattice and peak/valley neighborhoods select midplane points within a circle while the protein neighborhood will depend on the actual protein points and will be larger e.g. for the monomer
2. Perform the fit with fluid dimple center, and fixed decay to the midplane. This gives a smooth "dimple".
3. Measure the maximum curvature of the dimple over a set of points subject to one of the following conditions
	(1) Select the maximum curvature over the locations of the fitted points (neighborhood + buffer) even if the dimple center is elsewhere
	(2) Select the maximum curvature at the center of each dimple as long as the center is within one half-box vector of the neighborhood's center of mass (note that this is equivalent to "in the box" since I selected the neighborhood + buffer under PBCs)
	(3) Set the maximum curvature at the center of each dimple as long as the center is located in the (neighborhood + buffer) region
4. Filter the curvatures so that we only consider those with a |H_max| > 0.001 and |H_max| < 0.06
5. Report a H_max* according to one of the following rules
	(1) Unmodified H_max of all of the H_max values that satisfy the magnitude filtuer
	(2) Scale the H_max by the number of valid fits so that lower numbers of valid fits reduce the 
	
My preferred procedure is as follows:

--Use instantaneous peak/valley neighborhood
--Use a 15+ nm buffer
--Measure maximum curvature regardless of position
--Scale by valid frames, effectively giving an H_max* that implicitly describes the consistency of the surfaces
--Report both the peak and valley H_max* values which are balanced for the control and show much stronger curvature for the protein-membrane systems

Justification for my choices:

Using the peak/valley method is necessary to control for the number of points in the neighborhood (which varies with the size of the oligomer) and also make for a fair comparison to the control (i.e. you don't have to choose which protein points to import).
Larger buffers (above 10nm) give more consistent results across systems, and puts the absolute neighborhood+buffer size on par with the smaller buffers and the protein neighborhood method.
Measuring maximum curvature anywhere on the neighborhood+buffer means we are not restricted by a particular part of the dimple and instead we are just using it as a general smoothing tool. Incidentally, the protein systems will be reporting more dimple centers because of the focusing effect, so this choice actually helps distinguish the control better.
Scaling H_max by the number of fits puts an effective weight on the best fits, incorporating a measure of how focused the curvature is.
And finally, reporting the balance helps distinguish the control. Even though curvature has to sum to zero over the patch, the maximum curvature doesn't, so the imbalance in protein systems reflects the consolidation/focusing of positive curvature.

CHECK PBCS
DO RADAR PLOT 
and so on

end
'''

#---experiment parameters
params_expt = {
	'decay_z0' : 'zero',
	'cutoff' : 150,
	'fit_dynamics_type' : 'dynamic',
	'z_filter' : 'inclusive',
	'framecounts' : 500, 
	'geography' : '',
	}
	
#---sweeping parameters
params_sweeps = {
	'cutoff':[200,150,300],
	'geography':[
		'proteins',
		['control','v614-120000-220000-200'],
		['control','v612-75000-175000-200'],
		'max','min',
		'max_dynamic','min_dynamic',
		['lattice_grid',6],
		][3:7],
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
	'filter_type' : 'mod2',
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
			'cutoff' : 200,
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
			
#---PLOT PARAMETERS, META-SUMMARY
#-------------------------------------------------------------------------------------------------------------


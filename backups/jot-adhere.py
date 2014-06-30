#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')

#---MAIN
#-------------------------------------------------------------------------------------------------------------

director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
basedir = '/home/rpb/worker/membrain/'
struct_protein = 'adhesion-tests/protein-cg.gro'
struct_membrane = 'adhesion-tests/membrane-start-no-solvent.gro'
place_position = 'center'
z_sweep_start = -1.
z_sweep_end = 10.
z_sweep_interval = 1.

if 0:
	#---Load the membrane
	mset.load_trajectory((basedir+struct_membrane,basedir+struct_membrane),resolution='cgmd')
	mset.identify_monolayers(director_cgmd,startframeno=0)
	mset.calculate_midplane('name PO4',0,rounder=20.,interp='best')
	topxyz = array([mean(mset.universe.residues[i].selectAtoms(selector_cgmd).coordinates(),axis=0) for i in mset.monolayer_residues[0]])
	

#---Load the protein
if 0:
	prot = Universe(basedir+struct_protein,basedir+struct_protein)
	protpos = prot.selectAtoms('all').coordinates()
	grid_positions = []
	if place_position == 'center':
		place_xyz = [mset.universe.dimensions[0]/2.,mset.universe.dimensions[1]/2.,mset.universe.dimensions[2]]
		protpos = array([p+(place_xyz-mean(protpos,axis=0)) for p in protpos])
		grid_positions.append(place_xyz)

if 0:
	meshplot(topxyz,vecs=mset.vecs[0])
	#meshpoints(protpos,scale_factor=10,color=(1,1,1))
	dists = [norm(mean(protpos,axis=0)[0:2]-t[0:2]) for t in topxyz]
	closelipid = dists.index(min(dists))
	#meshpoints(topxyz[closelipid],scale_factor=20,color=(0,0,0))
	protein_configs = []
	for z_shift in arange(z_sweep_start,z_sweep_end,z_sweep_interval):
		shifter = [0,0,-1*(mean(protpos,axis=0)[2]-topxyz[closelipid][2])+10*z_shift]
		protein_configs.append([p+shifter for p in protpos])

if 0:
	meshplot(topxyz,vecs=mset.vecs[0])
	scale_factor=20
	color=(0,0,0)
	translate=[0.0,0.0,0.0]
	scale_factor=10
	color=(0,0,0)
	opacity=1
	resolution=8
	prot = protein_configs[0]
	data = prot
	delay = 0.5
	fig = mlab.gcf()
	#---meshpoint function
	if len(data) == 1 or (type(data) == numpy.ndarray and len(data) == 3):
		data = [data]
	if type(data) == list:
		data = array(data)
	X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
	if scale_factor == None:
		scale_factor = 3.
		scale_mode = 'none'
		s = Z
	elif type(scale_factor) == int:
		scale_factor = [scale_factor for i in range(len(data))]
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	else:
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	if shape(data)[1] != 3:
		data = unzipgrid(data,vecs=vecs)
	pts = mlab.points3d(X,Y,Z,s,scale_mode=scale_mode,mode='sphere',scale_factor=0.5,color=color,opacity=opacity,resolution=resolution)
	#---end
	for prot in protein_configs:
		data = prot
		#---meshpoint function
		if len(data) == 1 or (type(data) == numpy.ndarray and len(data) == 3):
			data = [data]
		if type(data) == list:
			data = array(data)
		X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
		if scale_factor == None:
			scale_factor = 3.
			scale_mode = 'none'
			s = Z
		elif type(scale_factor) == int:
			scale_factor = [scale_factor for i in range(len(data))]
			scale_mode = 'scalar'
			s = scale_factor
			scale_factor = 0.0005
		else:
			scale_mode = 'scalar'
			s = scale_factor
			scale_factor = 0.0005
		if shape(data)[1] != 3:
			data = unzipgrid(data,vecs=vecs)
		#---end
		pts.mlab_source.x = X
		pts.mlab_source.y = Y
		pts.mlab_source.z = Z
		pts.mlab_source.s = s
		raw_input('Enter to animate, cave-man style.')
		
if 1:
	for out_coords in protein_configs:
		fpin = open(basedir+struct_protein,'r')
		struct_protein_file = []
		for line in fpin:
			struct_protein_file.append(line)
		fpin.close()
		fpout = open(basedir+'/outfile.gro','w')
		fpout.write(struct_protein_file[0])
		fpout.write(struct_protein_file[1])
		for l in range(2,len(struct_protein_file)-1):
			fpout.write(struct_protein_file[l][0:21]+'  '+str('%.3f'%out_coords[l-2][0])+'  '+
				str('%.3f'%out_coords[l-2][1])+'  '+str('%.3f'%out_coords[l-2][0])+'\n')
		fpout.close()

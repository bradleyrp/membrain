#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')

#---PARAMETERS
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
Nxy = [2,2]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Load the membrane
mset.load_trajectory((basedir+struct_membrane,basedir+struct_membrane),resolution='cgmd')
mset.identify_monolayers(director_cgmd,startframeno=0)
topxyz = array([mean(mset.universe.residues[i].selectAtoms(selector_cgmd).coordinates(),axis=0) 
	for i in mset.monolayer_residues[0]])
vecs = mset.vec(0)

#---Define a simple lattice
lattice = [[vecs[0]*(i+1)/(Nxy[0]+1),vecs[1]*(j+1)/(Nxy[1]+1)] 
	for i in range(0,Nxy[1]) for j in range(0,Nxy[0])]

#---Load the protein
prot = Universe(basedir+struct_protein,basedir+struct_protein)
protpos = prot.selectAtoms('all').coordinates()

#---Compute lattice grid offsets for moving the protein into correct XY position
grid_positions = []
protein_configs = []
for lat in lattice:
	xypos = [lat[0],lat[1],vecs[2]]
	protpos = array([p+(xypos-mean(protpos,axis=0)) for p in protpos])
	grid_positions.append(xypos)
	#---Find the height of the closest lipid in XY
	dists = [norm(mean(protpos,axis=0)[0:2]-t[0:2]) for t in topxyz]
	closelipid = dists.index(min(dists))
	protein_configs_sweep = []
	for z_shift in arange(z_sweep_start,z_sweep_end,z_sweep_interval):
		shifter = [0,0,-1*(mean(protpos,axis=0)[2]-topxyz[closelipid][2])+10*z_shift]
		protein_configs_sweep.append([p+shifter for p in protpos])
	protein_configs.append(protein_configs_sweep)

#---Write the resulting frames
fplog = open('logfile','w')
for gp in range(len(grid_positions)):
	gridpos = grid_positions[gp]
	fplog.write('Position: '+str(gridpos[0:2])+'\n')
	for pp in range(len(protein_configs[gp])):
		out_coords = protein_configs[gp][pp]
		fpin = open(basedir+struct_protein,'r')
		struct_protein_file = []
		for line in fpin:
			struct_protein_file.append(line)
		fpin.close()
		fpout = open(basedir+'/protein_grid_'+str('%04d'%gp)+'_position_'+str('%04d'%pp)+'.gro','w')
		fpout.write(struct_protein_file[0])
		fpout.write(struct_protein_file[1])
		for l in range(2,len(struct_protein_file)-1):
			fpout.write(struct_protein_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
				str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.))+'\n')
		fpout.write(struct_protein_file[-1])
		fpout.close()
fplog.close()

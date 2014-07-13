#!/usr/bin/python -i

from membrainrunner import *
location = ''
execfile('locations.py')
execfile('plotter.py')

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

director_cgmd = ['name PO4','name C4A','name C4B']
selector_cgmd = 'name PO4'
basedir = '/home/rpb/worker/membrain/adhesion-tests/'
struct_protein = 'protein-cg.gro'
struct_membrane = 'membrane-start-no-solvent.gro'
top_membrane = 'system-bilayer.top'
top_protein = 'Protein_A.itp'
place_position = 'center'
z_sweep_start = -1.
z_sweep_end = 10.
z_sweep_interval = 1.
Nxy = [2,2]

#---method
write_combined_gro = True
write_separate_gro = False

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---load the membrane
mset.load_trajectory((basedir+struct_membrane,basedir+struct_membrane),resolution='cgmd')
mset.identify_monolayers(director_cgmd,startframeno=0)
topxyz = array([mean(mset.universe.residues[i].selectAtoms(selector_cgmd).coordinates(),axis=0) 
	for i in mset.monolayer_residues[0]])
vecs = mset.vec(0)

#---define a simple lattice
lattice = [[vecs[0]*(i+1)/(Nxy[0]+1),vecs[1]*(j+1)/(Nxy[1]+1)] 
	for i in range(0,Nxy[1]) for j in range(0,Nxy[0])]

#---load the protein
prot = Universe(basedir+struct_protein,basedir+struct_protein)
protpos = prot.selectAtoms('all').coordinates()

#---compute lattice grid offsets for moving the protein into correct XY position
grid_positions = []
protein_configs = []
for lat in lattice:
	xypos = [lat[0],lat[1],vecs[2]]
	protpos = array([p+(xypos-mean(protpos,axis=0)) for p in protpos])
	grid_positions.append(xypos)
	#---find the height of the closest lipid in XY
	dists = [scipy.linalg.norm(mean(protpos,axis=0)[0:2]-t[0:2]) for t in topxyz]
	closelipid = dists.index(min(dists))
	protein_configs_sweep = []
	for z_shift in arange(z_sweep_start,z_sweep_end,z_sweep_interval):
		shifter = [0,0,-1*(mean(protpos,axis=0)[2]-topxyz[closelipid][2])+10*z_shift]
		protein_configs_sweep.append([p+shifter for p in protpos])
	protein_configs.append(protein_configs_sweep)
	
#---load the bilayer gro file
fpin_bilayer = open(basedir+struct_membrane,'r')
struct_bilayer_file = []
for line in fpin_bilayer:
	struct_bilayer_file.append(line)
fpin_bilayer.close()

#---load the protein gro file
fpin = open(basedir+struct_protein,'r')
struct_protein_file = []
for line in fpin:
	struct_protein_file.append(line)
fpin.close()

#---write the resulting frames
if write_combined_gro:
	fplog = open('logfile','w')
	fpout = open(basedir+'/position-test.gro','w')
	for gp in range(len(grid_positions)):
		print 'writing test position '+str(gp)+' of '+str(len(grid_positions))
		gridpos = grid_positions[gp]
		fplog.write('Position: '+str(gridpos[0:2])+'\n')
		for pp in range(len(protein_configs[gp])):
			out_coords = protein_configs[gp][pp]
			fpout.write('bilayer sytem\n')
			fpout.write(str(len(struct_protein_file)+len(struct_bilayer_file)-6)+'\n')
			for l in range(2,len(struct_protein_file)-1):
				fpout.write(struct_protein_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
					str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.))+'\n')
			for l in range(2,len(struct_bilayer_file)-1):
				fpout.write(struct_bilayer_file[l])
			fpout.write(struct_bilayer_file[-1])
	fpout.close()
	fplog.close()
	#---write a reference structure
	fpout = open(basedir+'/position-test.start.gro','w')
	fpout.write('bilayer sytem\n')
	fpout.write(str(len(struct_protein_file)+len(struct_bilayer_file)-6)+'\n')
	for l in range(2,len(struct_protein_file)-1):
		fpout.write(struct_protein_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
			str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.))+'\n')
	for l in range(2,len(struct_bilayer_file)-1):
		fpout.write(struct_bilayer_file[l])
	fpout.write(struct_bilayer_file[-1])
	fpout.close()
elif write_separate_gro:
	fplog = open('logfile','w')
	for gp in range(len(grid_positions)):
		print 'writing test position '+str(gp)+' of '+str(len(grid_positions))
		gridpos = grid_positions[gp]
		fplog.write('Position: '+str(gridpos[0:2])+'\n')
		for pp in range(len(protein_configs[gp])):
			out_coords = protein_configs[gp][pp]
			fpout = open(basedir+'/protein_grid_'+str('%04d'%gp)+'_position_'+str('%04d'%pp)+'.gro','w')
			fpout.write('bilayer sytem\n')
			fpout.write(str(len(struct_protein_file)+len(struct_bilayer_file)-6)+'\n')
			for l in range(2,len(struct_protein_file)-1):
				fpout.write(struct_protein_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
					str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.))+'\n')
			for l in range(2,len(struct_bilayer_file)-1):
				fpout.write(struct_bilayer_file[l])
			fpout.write(struct_bilayer_file[-1])
			fpout.close()
	fplog.close()

#---construct a topology
fpin = open(basedir+top_membrane,'r')
basetop = []
for line in fpin:
	basetop.append(line)
fpin.close()



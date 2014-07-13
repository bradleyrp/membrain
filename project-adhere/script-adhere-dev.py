#!/usr/bin/python

import subprocess

'''

#---load the bilayer
fpin = open(basedir+struct_protein,'r')
struct_bilayer_file = []
for line in fpin:
	struct_bilayer_file.append(line)
fpin.close()


fpout = open(basedir+'/protein_grid_'+str('%04d'%gp)+'_position_'+str('%04d'%pp)+'.gro','w')
fpout.write(struct_protein_file[0])
fpout.write(struct_protein_file[1])
for l in range(2,len(struct_protein_file)-1):
	fpout.write(struct_protein_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
		str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.))+'\n')
fpout.write(struct_protein_file[-1])
fpout.close()

#---print the concatenated protein and bilayer system
natoms = len(struct_protein_file)+len(struct_bilayer_file)-6
subprocess.call(['echo','bilayer system'],stdout=open(basedir+'prep-combined.gro','w'))
subprocess.call(['echo',str(natoms)],stdout=open(basedir+'prep-combined.gro','w'))


echo $NATOMS >> prep-combined.gro
sed -e '1d' -e '2d' -e '$d' prep-protein-moved.gro >> prep-combined.gro
sed -e '1d' -e '2d' prep-membrane-start.gro >> prep-combined.gro
	
'''

#---construct a topology
fpin = open(basedir+top_membrane,'r')
basetop = []
for line in fpin:
	basetop.append(line)
fpin.close()

basetopnosol = [i for i in basetop 
	if (i[0:2] != 'W ' and
	i[0:4] != 'NA+ ' and
	i[0:4] != 'CL- ')]

fpout = open(basedir+'system.top','w')
for l in range(len(basetopnosol)):
	line = basetopnosol[l]
	if l == max([i for i in range(len(basetopnosol)) 
		if basetopnosol[i][0:8] == '#include']):
		fpout.write('#include "Protein_A.itp"\n')
	fpout.write(line)
fpout.write('Protein_A 1')
fpout.close()

# echo -e "keep 1\nname 0 PROTEIN\nr DOPC | r DOPS\nname 1 LIPIDS\nq\n" | make_ndx -f position-test.start.gro -o index-system.ndx
# grompp -c position-test.start.gro -o position-test.tpr -f input-md-in.mdp -p system.top -n index-system.ndx -maxwarn 10
# mdrun -s position-test.tpr -rerun position-test.xtc

#!/usr/bin/python -i

from membrainrunner import *
import numpy as np

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---Analysis parameters
location = 'light'
execfile('locations.py')
framecount = 30
skip = None

#---Selections
#---Note: this code also uses the following globals: sel_aamd_surfacer
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---Analysis plan
analysis_descriptors = [
	(['membrane-v509'],'all',director_symmetric,-1)]

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def transform(sel1,sel2):
	#---Following function translates and rotates three points using SVD.
	sel1b = sel1-mean(sel1,axis=0)
	sel2b = sel2-mean(sel2,axis=0)	
	A=dot(sel1b.T,sel2b)
	svdans = linalg.svd(A)
	m=np.diag((1,1,linalg.det(dot(svdans[2].T,svdans[0].T))))
	tmat=dot(dot((svdans[2].T),m),svdans[0].T)
	sel4 = array([list(dot(tmat,i)) for i in sel1b])
	#meshpoints(transpose(sel1b),opacity=0.2,color=(1,0,0))
	#meshpoints(transpose(sel2b),opacity=0.2,color=(0,1,0))
	#meshpoints(transpose(sel4),color=(1,1,1),opacity=1)
	return [mean(sel1-sel2,axis=0),tmat]
		
#---MAIN
#-------------------------------------------------------------------------------------------------------------

snapshotfile='/home/rpb/worker/worker-big/membrane-repository/look-v509/md.part0033.gro'
snapshotfile='/home/rpb/worker/worker-big/membrane-repository/look-v509/md.part0033.gro'
snapshotfile='/home/rpb/worker/repo-membrane/membrane-v531/md.part0003.8854.gro'
#snapshotfile='/home/rpb/worker-big/membrane-repository/membrane-v531/md.part0003.8854.gro'
#snapshotfile='/home/rpb/worker-big/membrane-repository/membrane-v532/md.part0003.9926.05.gro'
snapshotfile='/home/rpb/worker-big/membrane-repository/membrane-v532/md.part0003.9926.gro'

#---Load
if 1:
	ad = analysis_descriptors[0]
	starttime = time.time()
	print 'Starting analysis job.'
	for ad in analysis_descriptors:
		#---Load global variables for analysis function
		(tests,selector,director,trajno) = ad
		traj = trajectories[systems.index(tests[0])][trajno]
		#---Load the trajectory
		gro = structures[systems.index(tests[0])]
		basename = traj.split('/')[-1][:-4]
		sel_surfacer = sel_aamd_surfacer
		print 'Accessing '+basename+'.'
		mset.load_trajectory((snapshotfile,snapshotfile),resolution='aamd')
		director = ['name P','name C218','name C318']
		#mset.identify_monolayers(director)

#---Demonstrate the fitting
if 1:
	resnum_template = 10
	resnum_target = 20
	#---Load the template
	lipids = mset.universe.selectAtoms('resname PI2P')
	template_base = lipids.residues[resnum_template].selectAtoms('name C13 or name C14 or name C15').coordinates()
	template_location = np.mean(template_base,axis=0)
	template_base = template_base-template_location
	template_all = lipids.residues[resnum_template].selectAtoms('name C13 or name C14 or name C15 or name O4 or name P4 or name OP42 or name OP43 or name OP44').coordinates()
	template_all = template_all-template_location
	#---Load the target
	target_lipid = lipids.residues[resnum_target].selectAtoms('all').coordinates()
	target_base = lipids.residues[resnum_target].selectAtoms('name C12 or name C13 or name C14').coordinates()
	ring = lipids.residues[resnum_target].selectAtoms('name C13 or name C14 or name C15 or name C16 or name C12 or name C11').coordinates()
	moved = transform(target_base,template_base)
	overlay = dot(moved[1].T,template_all.T).T+moved[0]
	#---Demonstrate
	if 0:
		meshpoints(overlay,color=(1,1,1),opacity=0.5,scale_factor=4)
		meshpoints(concatenate((ring,lipids.residues[resnum_target].selectAtoms('name C1 or name C2 or name P').coordinates())),color=(0,1,0),scale_factor=3,opacity=0.5)

#---Testing
if 1:
	#---Find new phosphate group coordinates
	basecount = 3
	resnum_template = 10
	resnum_target = 20
	template_names = ['C13','C14','C15','H4','O4','P4','OP42','OP43','OP44']
	target_names = ['C12','C13','C14','H3','O3']
	#---Load the template
	lipids = mset.universe.selectAtoms('resname PI2P')
	template_base = []
	for name in template_names[0:5]:
		template_base.append(lipids.residues[resnum_template].selectAtoms('name '+name).coordinates()[0])
	template_base = array(template_base)
	template_location = np.mean(template_base,axis=0)
	template_base = template_base-template_location
	template_all = []
	for name in template_names:
		template_all.append(lipids.residues[resnum_template].selectAtoms('name '+name).coordinates()[0])
	template_all = array(template_all)-template_location
	#---Loop over target lipids and generate new coordinates
	new_coords0 = []
	for resnum in range(len(lipids.resnums())):
		target_lipid = lipids.residues[resnum].selectAtoms('all').coordinates()
		#target_base = lipids.residues[resnum].selectAtoms('name C12 or name C13 or name C14').coordinates()
		target_base = []
		for name in target_names[0:5]:
			target_base.append(lipids.residues[resnum].selectAtoms('name '+name).coordinates()[0])
		target_base = array(target_base)		
		moved = transform(target_base,template_base)
		overlay = dot(moved[1].T,template_all.T).T+moved[0]
		new_coords0.append(overlay[basecount:])

#---Generate new coordinates
if 1:
	#---Find new hydroxyl group coordinates
	basecount = 3
	resnum_template = 10
	resnum_target = 20
	template_names = ['C12','C13','C14','H3','O3','HO3']
	target_names = ['C13','C14','C15','H4','O4']
	#---Load the template
	lipids = mset.universe.selectAtoms('resname PI2P')
	template_base = []
	for name in template_names[0:5]:
		template_base.append(lipids.residues[resnum_template].selectAtoms('name '+name).coordinates()[0])
	template_base = array(template_base)
	template_location = np.mean(template_base,axis=0)
	template_base = template_base-template_location
	template_all = []
	for name in template_names:
		template_all.append(lipids.residues[resnum_template].selectAtoms('name '+name).coordinates()[0])
	template_all = array(template_all)-template_location
	#---Loop over target lipids and generate new coordinates
	new_coords1 = []
	for resnum in range(len(lipids.resnums())):
		target_lipid = lipids.residues[resnum].selectAtoms('all').coordinates()
		#target_base = lipids.residues[resnum].selectAtoms('name C13 or name C14 or name C15').coordinates()
		target_base = []
		for name in target_names[0:5]:
			target_base.append(lipids.residues[resnum].selectAtoms('name '+name).coordinates()[0])
		target_base = array(target_base)
		moved = transform(target_base,template_base)
		overlay = dot(moved[1].T,template_all.T).T+moved[0]
		new_coords1.append(overlay[basecount:])
	
#---Make the correct substitutions
if 1:
	change1_remove_name = ['H4','O4','P4','OP42','OP43','OP44']
	change2_remove_name = ['H3','O3','HO3']
	changes_resname = 'PI2P'
	change3_new_names = ['H3','O3','P3','OP32','OP33','OP34']
	change4_new_names = ['H4','O4','HO4']
	change3_precede = 'C13'
	#change3_precede_follow = 'HO2'
	change4_precede = 'C14'
	#change4_precede_follow = 'HO3'
	original = [line.strip('\n') for line in open(snapshotfile)]
	#---Remove old entries
	step1 = []
	for line in original:
		if line[5:10].strip() != changes_resname:
			step1.append(line)
		elif ((line[5:10].strip() == changes_resname) and
			(line[10:15].strip() not in change1_remove_name) and
			(line[10:15].strip() not in change2_remove_name)):
			step1.append(line)
	#---Add new entries
	step2 = []
	thislipid = 0
	parts = 0
	for line in step1:
		if ((line[5:10].strip() == changes_resname) and
			(line[10:15].strip() == change3_precede)):
			# moved this above to fix order problem
			step2.append(line+'\n')
			for j in range(len(new_coords0[thislipid])):
				step2.append(line[0:5]+changes_resname+change3_new_names[j].rjust(6)+'00000'+
					('%.3f'%(new_coords0[thislipid][j][0]/10.)).rjust(8)+''+
					('%.3f'%(new_coords0[thislipid][j][1]/10.)).rjust(8)+''+
					('%.3f'%(new_coords0[thislipid][j][2]/10.)).rjust(8)+'\n')
			print thislipid
			parts = parts + 1
			if parts == 2:
				thislipid = thislipid + 1
				parts = 0
		elif ((line[5:10].strip() == changes_resname) and
			(line[10:15].strip() == change4_precede)):
			# moved this above to fix order problem
			step2.append(line+'\n')
			for j in range(len(new_coords1[thislipid])):
				step2.append(line[0:5]+changes_resname+change4_new_names[j].rjust(6)+'00000'+
					('%.3f'%(new_coords1[thislipid][j][0]/10.)).rjust(8)+''+
					('%.3f'%(new_coords1[thislipid][j][1]/10.)).rjust(8)+''+
					('%.3f'%(new_coords1[thislipid][j][2]/10.)).rjust(8)+'\n')
			print thislipid
			parts = parts + 1
			if parts == 2:
				thislipid = thislipid + 1
				parts = 0
		else:
			step2.append(line+'\n')
	#---Write the structure
	fp = open('test.gro','w')
	for line in step2:
		fp.write(line)
	fp.close()	
	

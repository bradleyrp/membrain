#!/usr/bin/env python

#---MEMBRANEDATA DEFINITIONS
#-------------------------------------------------------------------------------------------------------------

'''
These definitions provide a flexible way to define new results types for membraindata.
The membraindata class allows for easy storage and lookup of different slices of heterogeneous analysis data.
'''

'''
Working list of non-deprecated pickledumps. 
Sort through these and add to the master list.
Generated these lists via "grep pickledump *"

List of actively-dumped pickles is:

script-aamd-apl.py
script-aamd-gr2d.py
script-aamd-ice-skate.py
script-aamd-structure.py
script-cgmd-structure.py
script-dimple-adv.py
script-dimple.py
script-mesoproc.py
script-topography.py

excluding

script-aamd-lipid-ion-voronoi.py
script-beta-double-layer.py

uncertain

script-coupling-batch.py
script-lipid-gr.py
script-lipid-micro.py
script-lipid-tilt.py

postprocs, uncertain

script-postproc-dimple.py
script-postproc-stressmaps.py
script-postproc-stress.py
script-postproc-tilefilter.py
script-postproc-tilefilter.py

I'm entering the still-in-use pkl types according to the most recent pkl which is given by
ls ../repo-pickles/pkl.* -ltrh

no MembraneData type in pkl.bilayer-coupling-sweep
headspan-headangle to pkl.headspan-headangle2.membrane-v515.a3-structure.s2-sim-compbio-md.part0007.20000-30000-10.pkl

'''

#---master dictionary which provides specifications for different data types
master_datatypes = {
	'gr2d':{
		'struct':{'frame':0,'monolayer':1,'lipid':2},
		'struct_opts':None,
		'description':'Basic object which actually contains 3D radial distribution function (i.e. g(r)) data for AAMD bilayers. The root data is a normalized histogram or probability curve of the likelihood of reaching another particle. Particle specs are in the pkl name and the notes.',
		},
	}
	
master_datatypes_deprecated = {
	'grvoronoi':{
		'struct':{'frame':0,'monolayer':1,'direction':2,'lipid':3},
		'struct_opts':{'direction' : {'z':0,'2d':1}},
		'description':'',
		},
	'topography3':{
		'struct':{'frame':0,'heights':1},
		'struct_opts':None,
		'description':'type: updated object which holds height distributions near certain neighborhoods (topography)',
		},
	'dimple3':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}},
		'description':'type: updated/advanced dimple fitting (dimple3 pkl objects)',
		},
	'dimple2':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}},
		'description':'type: updated/advanced dimple fitting (dimple2 pkl objects)',
		},
	'spanangle2':{
		'struct':{'frame':0,'resid':1,'type':2},
		'struct_opts':{'type': {'headspan':0,'headangle':1,'time':2,'resid':3}},
		'description':'',
		},
	'spanangle':{
		'struct':{'frame':0,'resid':1,'headspan':2,'headangle':3},
		'struct_opts':None,
		'description':'',
		},
	'ionskate':{
		'struct':{'type':0,'zone':1,'ion':2,'deltat':3,'start_frame':4},
		'struct_opts':{'type': {'mastermsd_zones':0,'distsxy':1,'distsz':2}},
		'description':'type: ionskate for storing ice-skating ion MSD decompositions Nb the mastermsd_zones array is very jagged so easy lookups might not work Nb this was basically scrapped due to memory issues',
		},
	'topography_transform':{
		'struct':{'frame':0},
		'struct_opts':None,
		'description':'type: topography_transform holds an array of minimum distance vs surface point height Nb this is an alternative to the topographycorrelate class',
		},
	'c0map':{
		'struct':{'frame':0},
		'struct_opts':None,
		'description':'c0map holds data in stressmap pickles after integrating the voxel-wise stress tensors',
		},
	'xyz_assoc':{
		'struct':{'frame':0,'monolayer':1,'lipid':2},
		'struct_opts':None,
		'description':'',
		},
	'dimple':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type' : {'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}},
		'description':'',
		},
	'lipid_positions':{
		'struct':{'frame':0,'type':1,'monolayer':2,'lipid':3},
		'struct_opts':{'type':{'position':0,'mindist':1}},
		'description':'',
		},
	'protdists':{
		'struct':{'frame':0,'monolayer':2,'lipid':3},
		'struct_opts':None,
		'description':'',
		},
	'tilts':{
		'struct':{'frame':0,'tail':1,'monolayer':2,'lipid':3},
		'struct_opts':{'tail':{'a':0,'b':1}},
		'description':'',
		},
	'surfnorms':{
		'struct':{'frame':0,'monolayer':1,'type':2,'lipid':3},
		'struct_opts':{'type':{'normals','areas'}},
		'description':'stores surface normal vectors',
		},
	'tilefilter_area_v1':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type':{'positive':0,'negative':1}},
		'description':'',
		},
	'tilt_deprecated':{
		'struct':{'frame':0,'type':1,'monolayer':2,'lipid':3},
		'struct_opts':{'type':{'angle':0,'area':1}},
		'description':'',
		},
	'dimple_for':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type':{'params':0,'maxhs':1,'maxhxys':2,'target_zones':3,'frameno':4}},
		'description':'retired method, formerly named "dimple" so be careful not to confuse them',
		},
	'dimple_filter':{
		'struct':{'frame':0,'type':1},
		'struct_opts':{'type':{'points':0,'domain':1,'protptsplus':2,'domainw':3,'frameno':4}},
		'description':'retire the following data type, representing a verbatim port from old method',
		},
	'cells':{
		'struct':{'frame':0,'monolayer':1,'type':2},
		'struct_opts':{'type':{'points':0,'voronoi':1,'areas':2}},
		'description':'',
		},
	'triangles':{
		'struct':{'frame':0,'monolayer':1,'type':2},
		'struct_opts':{'type':{'points':0,'lines':1,'n_perfect_triangles':2}},
		'description':'',
		},
	'gr':{
		'struct':{'frame':0,'monolayer':1},
		'struct_opts':None,
		'description':'',
		},
	}
	


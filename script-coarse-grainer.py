#!/usr/bin/python -i

from membrainrunner import *

from scipy.linalg import norm

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = None
location = ''
execfile('locations.py')

#---selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---analysis plan, single entry selection only
analysis_descriptors = [
	('membrane-v509','resname PI2P',director_symmetric,-1)]
analysis = analysis_descriptors[0]

#---MAPPING
#-------------------------------------------------------------------------------------------------------------

aamap = [
	['P5','OP52','OP53','OP54'],
	['P4','OP42','OP43','OP44'],
	['C5'],
	['C3'],
	['C1'],
	['P1','O11','O12','O13'],
	['C12','H12','OI1','C21','OI2'],
	['C13','H13A','H13B','OS1','C31','OS2'],
	['C22','H2R','H2S','C23','H3R','H3S','C24','H4R','H4S'],
	['C25','H51','C26','H61','C27','H7R','H7S','C28','H81'],
	['C29','H91','C210','H10R','H10S','C211','H111','C212','H121'],
	['C213','H13R','H13S','C214','H141','C215','H151','C216','H16R','H16S'],	
	['C217','H17R','H17S','C218','H18R','H18S','C219','H19R','H19S','C220','H20R','H20S'],	
	['C32','H2X','H2Y','C33','H3X','H3Y','C34','H4X','H4Y','C35','H5X','H5Y','C36','H6X','H6Y'],
	['C37','H7X','H7Y','C38','H8X','H8Y','C39','H9X','H9Y','C310','H10X','H10Y'],
	['C311','H11X','H11Y','C312','H12X','H12Y','C313','H13X','H13Y','C314','H14X','H14Y'],
	['C315','H15X','H15Y','C316','H16X','H16Y','C317','H17X','H17Y','C318','H18X','H18Y']
	]
	
aamap = [
	['P5','OP52','OP53','OP54'],
	['P4','OP42','OP43','OP44'],
	['C5'],
	['C3'],
	['C1'],
	['P1','O11','O12','O13'],
	['C12','OI1','C21','OI2'],
	['C13','OS1','C31','OS2'],
	['C22','C23','C24'],
	['C25','C26','C27','C28'],
	['C29','C210','C211','C212'],
	['C213','C214','C215','C216'],	
	['C217','C218','C219','C220'],	
	['C32','C33','C34','C35','C36'],
	['C37','C38','C39','C310'],
	['C311','C312','C313','C314'],
	['C315','C316','C317','C318']
	]

cgmap = ['PO5','PO4','RO1','RO2','RO3','POG','GL1','GL2','C1A','C2A','C3A','C4A','C5A','C1B','C2B','C3B',
	'C4B']

def make_selection(beads):
	selstring = 'resname PIP2 and (name '+beads[0]
	for i in beads[1:]:
		selstring += ' or name '+i
	selstring += ')'
	return selstring

beadmap = []
for i in aamap:
	beadmap.append(make_selection(i))
	
bonds = [
	[1,3],
	[2,4],
	[5,6],
	[6,7],
	[7,8],
	[7,9],
	[9,10],
	[10,11],
	[11,12],
	[12,13],
	[8,14],
	[14,15],
	[15,16],
	[16,17],
	[3,4],
	[4,5],
	[5,3],
	[1,6],
	[2,6],
	[1,2]
	]
	
angles = [
	[1,3,4],
	[3,4,2],
	[6,5,4],
	[1,3,5],
	[2,4,5],
	[3,5,6],
	[7,8,14],
	[5,6,7],     
	[6,7,8],
	[6,7,9],
	[8,14,15],
	[7,9,10],
	[9,10,11],
	[10,11,12],
	[11,12,13],
	[14,15,16],
	[15,16,17]
	]
	
dihedrals = [
	[1,3,4,2],
	[1,3,5,4],
	[3,5,4,2],
	[3,4,5,6],
	[2,4,5,6],
	[4,3,5,6],
	[1,3,5,6],
	[3,5,6,7],
	[5,6,7,9],
	[6,7,9,10],
	[7,9,10,11],
	[9,10,11,12],
	[10,11,12,13],
	[4,5,6,7],
	[6,7,8,14],
	[14,15,16,17],
	[8,14,15,16]
	]

#---update names from previous "PIP2" to new model "PI2P"
updated_names = True	
map_old_index_new = {
	'C1':'C11','H1':'H1','O1':'O1','C2':'C12','H2':'H2','O2':'O2','HO2':'HO2','C3':'C13','H3':'H3','O3':'O3',
	'HO3':'HO3','C4':'C14','H4':'H4','O4':'O4','P4':'P4','OP42':'OP42','OP43':'OP43','OP44':'OP44','C5':'C15',
	'H5':'H5','O5':'O5','P5':'P5','OP52':'OP52','OP53':'OP53','OP54':'OP54','H52':'H52','C6':'C16','H6':'H6',
	'O6':'O6','HO6':'HO6','O13':'O11','P1':'P','O11':'O13','O12':'O14','C11':'C1','H11A':'HA','H11B':'HB',
	'H12':'HS','C12':'C2','H13A':'HX','C13':'C3','H13B':'HY','OS2':'O31','OS1':'O32','C31':'C31','C32':'C32',
	'H2X':'H2X','H2Y':'H2Y','OI1':'O21','OI2':'O22','C21':'C21','C22':'C22','H2R':'H2R','H2S':'H2S',
	'C33':'C33','H3X':'H3X','H3Y':'H3Y','C34':'C34','H4X':'H4X','H4Y':'H4Y','C35':'C35','H5X':'H5X',
	'H5Y':'H5Y','C36':'C36','H6X':'H6X','H6Y':'H6Y','C37':'C37','H7X':'H7X','H7Y':'H7Y','C38':'C38',
	'H8X':'H8X','H8Y':'H8Y','C39':'C39','H9X':'H9X','H9Y':'H9Y','C310':'C310','H10X':'H10X','H10Y':'H10Y',
	'C311':'C311','H11X':'H11X','H11Y':'H11Y','C312':'C312','H12X':'H12X','H12Y':'H12Y','C313':'C313',
	'H13X':'H13X','H13Y':'H13Y','C314':'C314','H14X':'H14X','H14Y':'H14Y','C315':'C315','H15X':'H15X',
	'H15Y':'H15Y','C316':'C316','H16X':'H16X','H16Y':'H16Y','C317':'C317','H17X':'H17X','H17Y':'H17Y',
	'C318':'C318','H18X':'H18X','H18Y':'H18Y','H18Z':'H18Z','C23':'C23','H3S':'H3S','H3R':'H3R','C24':'C24',
	'H4S':'H4S','H4R':'H4R','C25':'C25','H51':'H51','C26':'C26','H61':'H61','C27':'C27','H7R':'H7R',
	'H7S':'H7S','C28':'C28','H81':'H81','C29':'C29','H91':'H91','C210':'C210','H10S':'H10S','H10R':'H10R',
	'C211':'C211','H111':'H111','C212':'C212','H121':'H121','C213':'C213','H13S':'H13S','H13R':'H13R',
	'C214':'C214','H141':'H141','C215':'C215','H151':'H151','C216':'C216','H16R':'H16R','H16S':'H16S',
	'C217':'C217','H17R':'H17R','H17S':'H17S','C218':'C218','H18R':'H18R','H18S':'H18S','C219':'C219',
	'H19R':'H19R','H19S':'H19S','C220':'C220','H20R':'H20R','H20S':'H20S','H20T':'H20T'}
#---N.b. I generated the list above with the following print trick, followed by find/replace convert to dict
'''
mset2 = MembraneSet()
mset2.load_trajectory(('/home/rpb/worker/maptune-v5-change-triangle/membrane-v8-800.10/md.part0003.gro',
	'/home/rpb/worker/maptune-v5-change-triangle/membrane-v8-800.10/md.part0003.gro'),resolution='aamd')
mset2.universe.selectAtoms('resname PIP2').residues[0]
tmp = mset2.universe.selectAtoms('resname PIP2').residues[0]
[tmp.atoms[i].name for i in range(150)]
for a in [tmp.atoms[i].name for i in range(150)]:
	print '[\''+a+'\',\'\'],'
'''

#---MAIN
#-------------------------------------------------------------------------------------------------------------

'''
procedure
	read in a bunch of PIP2 coordinates
	take the average angle and dihedral values for each
	choose standard force values for each of the values
	write the output to a new ITP file
		filter as necessary, to match the number of angles/bonds/dihedrals standard for this molecule
		also allow options for constraints, necessary for the ring group	
'''

#---CLASSES
#---basic frame,residue objects, solely for organizing inputs to the CGMD-AAMD mapping procedure below
#---Nb. the following codeblock is modified from the original maptune code, written by rpb
class frame(object):
	def __init__(self,num):
		self.num = num
		self.res = []
	def add_res(self,residue):
		self.res.append(residue)		
class residue(object):
	def __init__(self,resid,resname):
		self.resid = resid
		self.resname = resname
		self.namelist = []
		self.poslist = []
	def add_atom(self,name,pos):
		self.poslist.append(pos)
		self.namelist.append(name)

#---load
(test,selector,director,trajno) = analysis
traj = trajectories[systems.index(test)][trajno]
mset = MembraneSet()
gro = structures[systems.index(test)]
basename = traj.split('/')[-1][:-4]
sel_surfacer = sel_aamd_surfacer
print 'Accessing '+basename+'.'
mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),resolution='aamd')
nfr = 100
fr = [ frame(i) for i in range(nfr+1) ]
i = 0
for snap in mset.universe.trajectory:
	print snap
	for r in mset.universe.selectAtoms(selector).residues:
		newres = residue(r.resids()[0],r.name)
		for a in r:
			newres.add_atom(str(a.name),list(a.pos))
		fr[i].add_res(newres)
	i += 1
	if i > nfr:
		break

#---move the structures into simple frame/residue objects
frcg = [ frame(i) for i in range(nfr+1) ]
for t in range(nfr+1):
	for resnr in range(len(fr[t].res)):
		if updated_names:
			beadgroups = [[fr[t].res[resnr].namelist.index(map_old_index_new[j]) for j in i] for i in aamap]
		else:
			beadgroups = [[fr[t].res[resnr].namelist.index(j) for j in i] for i in aamap]
		respos = [array([fr[t].res[resnr].poslist[i] for i in j]).mean(0) for j in beadgroups]
		newres = residue(r.resids()[0],'PIP2')
		for i in range(len(respos)):
			newres.add_atom(cgmap[i],list(respos[i]))
		frcg[t].add_res(newres)

#---extract parameters from simple frame/residue objects		
bondd = []
for t in range(nfr):
	framebondd = []
	for resnr in range(len(frcg[t].res)):
		framebondd.append(\
			[norm([a - b for a, b in zip(i[0], i[1])]) for i in \
			[[frcg[t].res[resnr].poslist[i-1] for i in j] for j in bonds]])
	bondd.append(framebondd)
angled = []
for t in range(nfr):
	frameangled = []
	for resnr in range(len(frcg[t].res)):
		#---atomized calculation, added to list below
		#i = [[frcg[0].res[0].poslist[i-1] for i in j] for j in angles][0]
		#[v1,v2] = transpose([[a-b,c-b] for a,b,c in zip(i[0], i[1],i[2])])
		#arccos(dot(v1,v2)/norm(v1)/norm(v2))
		frameangled.append(\
			[arccos(dot(v[0],v[1])/norm(v[0])/norm(v[1]))/pi*180 for v in \
			[transpose([[a-b,c-b] for a,b,c in zip(i[0], i[1],i[2])]) \
			for i in [[frcg[t].res[resnr].poslist[i-1] for i in j] for j in angles]]] \
			)
	angled.append(frameangled)
dihedd = []
for t in range(nfr):
	framedihedd = []
	for resnr in range(len(frcg[t].res)):
		snapresult = []
		for k in array([[frcg[t].res[resnr].poslist[i-1] for i in j] for j in dihedrals]):
			v10 = k[0]-k[1]
			v12 = k[2]-k[1]
			v23 = k[3]-k[2]
			v21 = k[1]-k[2]
			b1 = cross(v10,v12)
			b2 = cross(v23,v21)
			p1 = b1-dot(v12,b1)*v12/norm(v12)
			p2 = b2-dot(v12,b2)*v12/norm(v12)
			snapresult.append(arccos(dot(p1,p2)/norm(p1)/norm(p2))/pi*180)
		framedihedd.append(snapresult)
	dihedd.append(framedihedd)
	
#---compute structure parameters
bond_aa_means = list(array(bondd).mean(1).mean(0))
bond_aa_stds = list(array(bondd).std(0).std(0))
angle_aa_means = list(array(angled).mean(1).mean(0))
angle_aa_stds = list(array(angled).std(0).std(0))
dihed_aa_means = list(mean(mean(ma.masked_array(dihedd,isnan(dihedd)),axis=1),axis=0).filled(np.nan))
dihed_aa_stds =  list(std(std(ma.masked_array(dihedd,isnan(dihedd)),axis=1),axis=0).filled(np.nan))

#---filter these according to custom rules

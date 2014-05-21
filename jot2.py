
#--------	WORKING ON RERACKER

self = mset

atomdirectors = director
startframeno=0

'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
print 'status: identifying monolayers'
print 'status: moving to frame '+str(startframeno)
self.gotoframe(startframeno)
pointouts = self.universe.selectAtoms(atomdirectors[0])
pointins = [self.universe.selectAtoms(atomdirectors[j]).coordinates() 
	for j in range(1,len(atomdirectors))]
whichlayer = [0 if i > 0.0 else 1 for i in pointouts.coordinates()[:,2] - mean(pointins,axis=0)[:,2]]
monos = []
#---monos separates the lipids by absolute index into monolayers from index zero
monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 0])
monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 1])
#---monolayer rerack hack if some residue IDs are missing
#---Nb this may affect the tilter, mesher, identify_residues, and batch_gr functions so beware
if (max(monos[0]+monos[1])-min(monos[0]+monos[1])) != len(monos[0]+monos[1])-1:
	print 'warning: resorting the monolayer indices because there is a mismatch'
	#---reracker is a sorted list of all of the absolute indices
	reracker = list(sort(monos[0]+monos[1]))
	#---monos_rerack is a copy of reracker in relative indices which is separated into monolayers 
	monos_rerack = [[reracker.index(i) for i in monos[m]] for m in range(2)]
	self.monolayer_residues_abs = monos
	#---resids_reracker is in absolute units
	self.resids_reracker = reracker
	#---monolayer_residues is in relative units
	self.monolayer_residues = monos_rerack
else: self.monolayer_residues = monos
if len(monos[0]) != len(monos[1]):
	print 'warning: there is a difference in the number of lipids per monolayer'

print 'monos'
print [[min(i),max(i)] for i in monos]

print 'reracker'
print min(reracker),max(reracker)

print 'monosrerack'
print [[min(i),max(i)] for i in monos_rerack]

############ SELECTOR

selector = residues

'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
self.monolayer_by_resid = []
self.resnames = []
self.resids = []
print 'selector = '
print selector
#---selector holds the proposed resnames
for sel in selector:
	selection = self.universe.selectAtoms('resname '+sel)
	if len([i-1 for i in selection.resids()]):
		#---resids holds the absolute residue numbers
		self.resids.append([i-1 for i in selection.resids()])
		self.resnames.append(sel)
		
print 'self.resids'
print [[min(i),max(i)] for i in self.resids if len(i)>0]
		
#---when re-racker has been defined by identify_monolayers then provide distinct residue numberings
if self.resids_reracker != []: 
	self.resids_abs = self.resids
	newresids = [[self.resids_reracker.index(i) for i in j] for j in self.resids]
	self.resids = newresids
	#self.resids = [[flatten(self.monolayer_residues).index(i) for i in j] for j in self.resids]

else:
	self.resids_abs = self.resids
'''
#---absolute vs relative numbering is tricky here
if self.resids_reracker != []:
	#---note: errors observed here with APL script so rpb made the following change
	#---note: this is an indexing nightmare
	residue_ids = [[self.resids_reracker.index(i) for i in j] for j in self.resids]
	#residue_ids = [[self.monolayer_residues.index(i) for i in j] for j in self.resids]
else: residue_ids = self.resids
'''
#---monolayer_residues is in relative indices
for monolayer in self.monolayer_residues:
	monolayer_resids = []
	#for resids in residue_ids:
	#---self.resids is in relative indices
	for resids in self.resids:
		monolayer_resids.append(list(set.intersection(set(monolayer),set(resids))))
	#---monolayer_by_resid is in relative indices
	self.monolayer_by_resid.append(monolayer_resids)
if self.resids_reracker == []:
	self.monolayer_by_resid_abs = self.monolayer_by_resid
else:
	#---pack absolute residue indices into monolayer_by_resid_abs
	self.monolayer_by_resid_abs = \
		[[[self.resids_reracker[r] for r in restype] 
		for restype in mono] for mono in self.monolayer_by_resid]


#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')

'''
pickleset = ['pkl.lipid-ion-voronoi.v533.md.part0003.2000-5000-20.P35P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v534.md.part0003.2000-5000-20.P35P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v531.md.part0010.25500-27500-4.PI2P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v532.md.part0010.25500-27500-4.PI2P-to-ions.pkl']
names = ['PI(3,5)P2, Mg','PI(3,5)P2, Ca','PI(4,5)P2, Mg','PI(4,5)P2, Ca']
pairselects = ['P35P P-to-MG','P35P P-to-CA','PI2P P-to-MG','PI2P P-to-CA']
'''

analysis_descriptors = [
	('pkl.lipid-ion-voronoi.v533.md.part0003.2000-5000-20.P35P-to-ions.pkl','PI(3,5)P2, Mg','P35P P-to-MG',
	'P35P'),
	('pkl.lipid-ion-voronoi.v534.md.part0003.2000-5000-20.P35P-to-ions.pkl','PI(3,5)P2, Ca','P35P P-to-CA',
	'P35P'),
	('pkl.lipid-ion-voronoi.v531.md.part0010.25500-27500-4.PI2P-to-ions.pkl','PI(4,5)P2, Mg','PI2P P-to-MG',
	'PI2P'),
	('pkl.lipid-ion-voronoi.v532.md.part0010.25500-27500-4.PI2P-to-ions.pkl','PI(4,5)P2, Ca','PI2P P-to-CA',
	'PI2P')
	]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

msets = []
pickleset = [i[0] for i in analysis_descriptors]
for i in range(len(pickleset)):
	mset = unpickle(pickles+pickleset[i])
	msets.append(mset)

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlim((0,20))
for sys in range(len(msets)):
	mset = msets[sys]
	(name,pairselect,resname) = analysis_descriptors[sys][1:]
	md = mset.getdata(pairselect)
	dists2d = md.get(['monolayer',0,'lipid',resname,'direction','2d'],flat=False)
	dists1d = md.get(['monolayer',0,'lipid',resname,'direction','z'],flat=False)
	#dists1dflat = md.get(['monolayer',0,'lipid',resname,'direction','z'])
	absdists = [list(sqrt(array(dists2d[i][j])**2+array(dists1d[i][j])**2)) for i in range(len(dists1d)) for j in range(len(dists1d[i]))]
	absdists = [i for j in absdists for i in j]
	hist0 = numpy.histogram(absdists,bins=30,normed=True,range=(0,20))
	c = clrs[sys%len(clrs)]
	ax.plot(hist0[1][1:],hist0[0],'o-',lw=2,color=c,label=name)
plt.legend()
plt.show()
	

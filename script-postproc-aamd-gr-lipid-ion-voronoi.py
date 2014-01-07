#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')

pickleset = ['pkl.lipid-ion-voronoi.v533.md.part0003.2000-5000-20.P35P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v534.md.part0003.2000-5000-20.P35P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v531.md.part0010.25500-27500-4.PI2P-to-ions.pkl',
	'pkl.lipid-ion-voronoi.v532.md.part0010.25500-27500-4.PI2P-to-ions.pkl']
names = ['PI(3,5)P2, Mg','PI(3,5)P2, Ca','PI(4,5)P2, Mg','PI(4,5)P2, Ca']
pairselects = ['P35P P-to-MG','P35P P-to-CA','PI2P P-to-MG','PI2P P-to-CA']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

msets = []
for i in range(len(pickleset)):
	mset = unpickle(pickles+pickleset[i])
	msets.append(mset)

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlim((0,20))
for sys in range(len(msets)):
	mset = msets[sys]
	md = mset.getdata(pairselects[sys])
	dists2d = md.get(['monolayer',0,'lipid','P35P','direction','2d'],flat=False)
	dists1d = md.get(['monolayer',0,'lipid','P35P','direction','z'],flat=False)
	absdists = [list(sqrt(array(dists2d[i][j])**2+array(dists1d[i][j])**2)) for i in range(len(dists1d)) for j in range(len(dists1d[i]))]
	absdists = [i for j in absdists for i in j]
	dists1dflat = md.get(['monolayer',0,'lipid','P35P','direction','z'])
	hists = []
	hist0 = numpy.histogram(absdists,bins=30,normed=True,range=(0,20))
	hists.append(hist0)
	c = clrs[sys%len(clrs)]
	ax.plot(hist0[1][1:],hist0[0],'o-',lw=2,color=c,label=names[sys])
plt.legend()
plt.show()
	
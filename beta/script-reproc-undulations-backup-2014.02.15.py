#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

mset = unpickle(pickles+'pkl.structures.membrane-v614.s9-lonestar.120000-220000-200.pkl')

#---calculate 2D spectrum
mset.calculate_undulation_spectrum(removeavg=0,redundant=1)
m,n = shape(mset.uqraw)[1:]
spec_center = [int(i/2) for i in shape(mset.uqraw)[1:]]
lxy = array([mset.vec(i) for i in mset.surf_index])/mset.lenscale
xv,yv = meshgrid(range(m),range(n))
qs = [(sqrt((2*pi*(array(xv)-spec_center[0])/lxy[f][0])**2+
	(2*pi*(array(yv)-spec_center[1])/lxy[0][1])**2)) 
	for f in range(len(lxy))]
hqhq2d = mean(array(1.*(abs(array(mset.uqraw))/double(m*n))**2),axis=0)
qmag2d = mean(qs,axis=0)
plt.imshow(array(hqhq2d).T,interpolation='nearest',origin='lower');plt.show()

#---convert to 1D spectrum
spec1d = array([qmag2d,hqhq2d]).T.reshape(m*n,2)
qmagfilter = [10**-10,0.5]
specfilt = array(filter(lambda x: x[0] >= qmagfilter[0] and x[0] <= qmagfilter[1], spec1d))
lenscale = 10.
[bz,az]=numpy.polyfit(log(specfilt[:,0]),log(specfilt[:,1]),1)

#---compute
area = double(mean([i[0]*i[1] for i in mset.vecs])/lenscale**2)
print 'q-magnitude scaling = '+str(bz)
kappa = 1/exp(az)/area
print 'kappa = '+str(kappa)
specplot = array(filter(lambda x: x[0] >= 1./10.*qmagfilter[0] and x[0] <= qmagfilter[1]*10., spec1d))

#---plot
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(spec1d[:,0],spec1d[:,1],'o')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(specplot[:,0],[exp(az)*(i**bz) for i in specplot[:,0]])
plt.show()

'''
	#---fitting the best points
	specsort = spec1d[np.lexsort((spec1d[:,1],spec1d[:,0]))]
	best_rmsd = 10**10
	best_endpost = 0
	for endpost in range(0,len(specsort)/32):
		specfilter = specsort[endpost:len(specsort)/4]
		[bz,az] = numpy.polyfit(log(specfilter[:,0]),log(specfilter[:,1]),1)
		if sum([(exp(az)*(specfilter[i,0]**bz)-specfilter[i,1])**2 
			for i in range(len(specfilter))])/len(specfilter[:,0]) < best_rmsd:
			best_rmsd = sum([(exp(az)*(specfilter[i,0]**bz)-specfilter[i,1])**2 
				for i in range(len(specfilter))])/len(specfilter[:,0])
			best_endpost = endpost
	qmagfilter = [specsort[best_endpost,0],specsort[len(specsort)/32],0]]

'''

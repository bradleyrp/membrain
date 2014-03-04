#!/usr/bin/python

#---look at the MSD of a 1D brownian particle with a harmonic restoring force
#---demonstrates that when MSD curves return to zero at high durations they represent something that is effectively bound

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

if 1:
	lims = [-1,1]
	maxt = 500
	nparts = 10
	bias = True
	parts = []
	for i in range(nparts):
		xs = [0]
		x = xs[0]
		for step in range(1,maxt):
			print 'step = '+str(step)+' x = '+str(x)
			if bias:
				restore = 1./50.*(x-xs[0])**2
				x += (lims[1]-lims[0])*np.random.random_sample()+lims[0] + (-1*restore if x > 0 else restore)
			else:
				x += (lims[1]-lims[0])*np.random.random_sample()+lims[0]
			xs.append(x)
		parts.append(xs)
	parts = np.array(parts)
	dists = []
	for p in range(len(parts)):
		print p
		part = parts[p]
		dists.append([[norm(part[s+d]-part[s])**2 for s in range(maxt-d)] for d in range(maxt)])
fig = plt.figure()
ax = plt.subplot(111)
for dist in dists:
	ax.plot([np.mean(i) for i in dist]);
ax.set_xscale('log')
ax.set_yscale('log')
plt.show()

	
'''
if 0:
	for part in parts:
		plt.plot(array(part)**2)
	plt.show()
dists = []
for i in range(nparts)):
	print i
	dists.append()
k = 10
dists = [[norm(parts[k,s:s+d]-parts[k,s])**2 for s in range(maxt-d)] for d in range(maxt)]
meandists = 
if 0:
	fig = plt.figure()
	ax = plt.subplot(111)
	meanxs = np.mean(parts,axis=0)
	ax.plot(array(meanxs))
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.show()
'''

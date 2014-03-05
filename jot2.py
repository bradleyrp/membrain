#!/usr/bin/python

#for rawdat in tilefiltdat:
#	plt.plot([float(sum(rawdat[rawdat[:,0]<i,1]>0.))/sum(rawdat[:,0]<i) for i in range(10,300,10)])
#plt.show()
if 0:
	areacurves = []
	for rawdat in tilefiltdat:
		areacurves.append([float(sum(rawdat[rawdat[:,0]<i,1]>0.))/sum(rawdat[:,0]<i) for i in range(10,300,10)])
	areacurves = array(areacurves)

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(mean(areacurves,axis=0))
ax.plot(mean(areacurves,axis=0)+std(areacurves,axis=0))
ax.plot(mean(areacurves,axis=0)-std(areacurves,axis=0))
plt.show()


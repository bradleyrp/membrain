#!/usr/bin/python -i

from membrainrunner import *
execfile('locations.py')

xdat,ydat = unpickle(pickles+'tmp.v616.late.sharp.spectrum.pkl')

ax = plt.subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim((0.1,10))
ax.set_xlim((0.01,10))			
ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
ax.grid(True)

xdat,ydat = unpickle(pickles+'tmp.v616.late.sharp.spectrum.pkl')
ax.plot(xdat,ydat,'r-',lw=2)
xdat,ydat = unpickle(pickles+'tmp.v616.late.dull.spectrum.pkl')
ax.plot(xdat,ydat,'b-',lw=2)

plt.show()

#!/usr/bin/python -i

from membrainrunner import *
execfile('plotter.py')

if 1:
	mv514=unpickle('pkl.gr-vornoi.membrane-v514.md.part0008.skip50.pkl')
	mv515=unpickle('pkl.gr-vornoi.membrane-v515.md.part0008.skip50.pkl')

histplots = plot_gr_voronoi_compare(mv514,mv515,label0=[mv514.avail()[0]],label1=[mv515.avail()[0]],colors=[0,1])
import numpy as np
from scipy.optimize import curve_fit

def func2(x,c1,c2,c3):
	return c1*exp(2.1*(exp(-c2*x)))

def func(x,c1,c2,c3):
	return c1*exp(c3*(exp(-c2*x)))
	
#blue
y = histplots[0][0]
x = histplots[0][1][:-1]
popt, pcov = curve_fit(func2, x, y)
print popt
plot(x,func2(x,popt[0],popt[1],popt[2]), color=clrs[1], linewidth=3)

y = histplots[1][0]
x = histplots[1][1][:-1]
popt, pcov = curve_fit(func, x, y)
print popt
plot(x,func(x,popt[0],popt[1],popt[2]),color=clrs[0], linewidth=3)
plt.show()

	

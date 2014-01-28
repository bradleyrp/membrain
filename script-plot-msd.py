#!/usr/bin/python -i

import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *

lines = []
fp = open('msd-v532.xvg','r')
for line in fp:
	if line[0] != '#' and line[0] != '@':
		lines.append([float(i) for i in line.split()])
fp.close()

lines = array(lines)
tdatall = log10(lines[1:])
tdat = log10(lines[int(len(tdatall)*.1)+2:int(len(tdatall)*.9)+1])
[bz,az] = polyfit(tdat[:,0],tdat[:,1],1)
print bz
fitted = [bz*i+az for i in tdatall[:,0]]
plt.plot(tdatall[:,0],tdatall[:,1],'-',c='b',lw=2)
plt.plot(tdat[:,0],tdat[:,1],'-',c='r',lw=2)
plt.plot(tdatall[:,0],fitted,'-',c='k',lw=1,alpha=0.5)
plt.show()

#!/usr/bin/python

import matplotlib as mpl
fullans = ms.fullans
n2,m2 = ms.m2,ms.n2
rowcombo = mean([array([[fullans[1][fr][i+(j-1)*m2] for j in range(n2)] for i in range(m2)])*fullans[0][fr]**2 for fr in range(m2*n2)],axis=0)
plt.imshow(abs(rowcombo).T,interpolation='nearest',origin='lower',norm=(mpl.colors.LogNorm() if 0 else None),cmap=mpl.cm.binary);plt.show()


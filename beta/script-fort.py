#!/usr/bin/python

import numpy
import curvcalc
dpts = numpy.array([[i for i in range(10)] for j in range(3)])
print dpts
print type(dpts)
print len(dpts[0])
curvcalc.curvcalc_solo(len(dpts),dpts)


import numpy
import memint
import time
import math

pts = numpy.array([1.12312,2,3,4,5,6,7,8,9,0,0,0])
#pts = numpy.array([[1,2,3],[4,5,6],[7,8,9],[0,0,0]])
#pts = numpy.zeros((3,2))
#ptstri = numpy.array([10,11,12,10,11,12,10,11,12])
memint.init(300,300)
memint.addvertex(0,numpy.array([1.12312,2.,3.]),numpy.array([6,7,8]),numpy.array([6,7]))
#memint.check(0)
#f = numpy.vectorize(memint.func)
#print f(pts)


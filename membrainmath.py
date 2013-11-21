#!/usr/bin/python

#---Mean curvature of a fluid-axis 2D Gaussian dimple

'''
# Calculation instructions follow

from sympy import *

z0 = Symbol('z0')
c0 = Symbol('c0')
x0 = Symbol('x0')
y0 = Symbol('y0')
sx = Symbol('sx')
sy = Symbol('sy')
th = Symbol('th')
x = Symbol('x')
y = Symbol('y')

fgauss2dfixed = z0+c0*exp(-(x-x0)**2/2/sx**2)*exp(-(y-y0)**2/2/sy**2)
fgauss2d = z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2/sx**2)*exp(-(-(x-x0)*sin(th)+(y-y0)*cos(th))**2/2/sy**2)

zx = diff(fgauss2d,x,1)
zy = diff(fgauss2d,y,1)
zxx = diff(zx,x,1)
zyy = diff(zy,y,1)
zxy = diff(zx,y,1)

meancurv = ((1+zx**2)*zyy+(1+zy**2)*zxx-2*zx*zy*zxy)/(1+zx**2+zy*2)*(3/2)

from sympy import *
'''

def gauss2dh(params,x,y):
	'''Mean curvature of a fluid-axis 2D Gaussian dimple.'''
	z0,c0,x0,y0,sx,sy,th = params
	return (-(2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - 2*c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)*(-c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2)*(c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sy**2 - c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sx**2 - c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/(sx**2*sy**2) + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sx**4) + ((c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)**2 + 1)*(-c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sy**2 + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sx**2 + 2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sx**4) + ((-c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2)**2 + 1)*(-c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sy**2 + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sx**2 - 2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sx**4))/(-2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - 2*c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2 + (c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)**2 + 1)
	

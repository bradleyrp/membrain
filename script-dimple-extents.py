#!/usr/bin/python

interact = True
from membrainrunner import *
execfile('locations.py')

from scipy import optimize

def gauss2d(params,x,y):
	'''Two-dimensional Gaussian height function with fluid axis of rotation.'''
	z0,c0,x0,y0,sx,sy,th = params
	#return a+b*exp(-(((x-c1)*cos(t1)+(y-c2)*sin(t1))/w1)**2-((-(x-c1)*sin(t1)+(y-c2)*cos(t1))/w2)**2)
	return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
		(y-y0)*cos(th))**2/2./sy**2)

def gauss2dh(params,x,y):
	z0,c0,x0,y0,sx,sy,th = params
	return ((2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - 2*c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)*(-c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2)*(c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sy**2 - c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sx**2 - c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/(sx**2*sy**2) + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/sx**4) + ((c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)**2 + 1)*(-c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sy**2 + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sx**2 + 2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sx**4) + ((-c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2)**2 + 1)*(-c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sy**2 + c0*((-x + x0)*sin(th) + (y - y0)*cos(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)**2/sy**4 - c0*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sx**2 - 2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)*cos(th)/(sx**2*sy**2) + c0*((x - x0)*cos(th) + (y - y0)*sin(th))**2*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)**2/sx**4))/(-2*c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sy**2 - 2*c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sx**2 + (c0*((-x + x0)*sin(th) + (y - y0)*cos(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*sin(th)/sy**2 - c0*((x - x0)*cos(th) + (y - y0)*sin(th))*exp(-((x - x0)*cos(th) + (y - y0)*sin(th))**2/(2*sx**2))*exp(-((-x + x0)*sin(th) + (y - y0)*cos(th))**2/(2*sy**2))*cos(th)/sx**2)**2 + 1)

fig = plt.figure(figsize=(8,6))
symbollist = ['o','s','>','<','^']	
xlim = 20
sy = 2
ax = plt.subplot(121)
sigxs = [2,4,8,16]
for sxi in range(len(sigxs)):
	sx = sigxs[sxi]
	params = [0,0.04,0,0,sx,2,0]
	xpts = linspace(-xlim,xlim,100)
	hxpts = [-1*gauss2dh(params,x,0) for x in xpts]
	hypts = [-1*gauss2dh(params,0,y) for y in xpts]
	ax.plot(xpts,hxpts,'-',c='r',lw=2,alpha=1-float(sxi)/len(sigxs),label='long, '+str(sx))
	ax.plot(xpts,hypts,'-',c='b',lw=2,alpha=1-float(sxi)/len(sigxs),label='long, '+str(sx))
ax.set_xlim((-xlim,xlim))
ax.set_ylim((-0.008,0.025))
ax.set_xticks([-10,0,10])
ax.axhline(y=0,xmin=0,xmax=1,lw=2)
ax.grid(True)
ax.set_title(r'dimple curvature',fontsize=fsaxtitle)
plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
ax.set_ylabel(r'$H\,\mathrm{({nm}^{-1})}$',fontsize=fsaxlabel)
ax.set_xlabel(r'$x\,\mathrm{(nm)}$',fontsize=fsaxlabel)
sigys = linspace(2,4,3)
sigys = arange(2,4+0.5,0.25)
isotrop_pts = []
for syi in range(len(sigys)):
	sy = sigys[syi]
	ax = plt.subplot(122)
	collectx = []
	collecty = []
	sigmas = [2,4,8,16]
	sigmas = linspace(sy,6,10)
	init = 0
	for s in range(len(sigmas)):
		sx = sigmas[s]
		if s == 0: initx = sx/2
		if s == 0: inity = sy/2
		params = [0,0.1,0,0,sx,sy,0]
		testx = lambda x : gauss2dh(params,x,0)**2
		res = optimize.fmin(testx,initx,ftol=0.0001)
		initx = res[0]
		collectx.append(res[0])
		testy = lambda y : gauss2dh(params,0,y)**2
		res = optimize.fmin(testy,inity,ftol=0.0001)
		init = res[0]
		collecty.append(res[0])
	isotrop_pts.append(collectx[0])
	if sy in sigys[::4]:
		syi2 = where(sigys[::4]==sy)[0][0]
		ax.plot(sigmas,collecty,symbollist[syi2]+'-',c='b')
		ax.plot(sigmas,collectx,symbollist[syi2]+'-',c='r')
ax.plot(sigys,isotrop_pts,'-',lw=2)
ax.grid(True)
ax.set_xlim((1.5,6.5))
ax.set_ylim((1,20))
ax.set_xticks([2,4,6])
ax.set_title(r'curvature extents',fontsize=fsaxtitle)
plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
ax.set_xlabel(r'$\sigma_x\,\mathrm{(nm)}$',fontsize=fsaxlabel)
ax.set_ylabel(r'$(\sigma_{x}^{*},\sigma_{y}^{*})\,\mathrm{(nm)}$',fontsize=fsaxlabel)
plt.tight_layout()
plt.savefig(pickles+'fig-dimple-extents-compare.png',dpi=300,bbox_inches='tight')
plt.show()


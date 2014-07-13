#!/usr/bin/python

from membrainrunner import *
execfile('locations.py')

import sys,os,re,time
import datetime
import psycopg2
import psycopg2.extras

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---settings
plotshow = False

#---INCLUDES
#-------------------------------------------------------------------------------------------------------------

#---connect
if 'conn' not in globals(): 
	try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	except: raise Exception('except: cannot connect to database')
	cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

if 1:
	#---retrieve data
	cur.execute('SELECT kappa_apparent,mesosims.kappa,mesosims.c_0,mesosims.callsign '+\
		'FROM dataref_structure,mesosims '+\
		'WHERE mesosims.id=dataref_structure.parent_mesosims')
	dat = cur.fetchall()
	kappas_app,kappas,curvs = array([i[:3] for i in dat]).T
	callsigns = [i[3] for i in dat]
	colorset = [clrs[list(set(callsigns)).index(i)] for i in callsigns]
	#---compute means and errors over the unique curvature values
	curvs_sorted = sort(list(set(curvs)))
	inds = [list(where(curvs==c)[0]) for c in curvs_sorted]
	means = [mean([(kappas_app-kappas)[i] for i in ilist]) for ilist in inds]
	stddevs = [std([(kappas_app-kappas)[i] for i in ilist]) for ilist in inds]
	#---isolate zero curvatures
	inds_zero = where(curvs==0)[0]
	#---deprecated plot of kappa versus kappa_apparent
	if 0:
		ax = plt.subplot(121)
		ax.scatter(kappas_app,kappas,c=colorset,s=30,lw=0)
		ax.scatter(kappas_app[inds_zero],kappas[inds_zero],c='k',s=40,lw=0)
		#ax.plot(range(40),range(40),'k-',lw=2)
		ax.grid(True)
	fig = plt.figure(figsize=(6,6))
	ax = plt.subplot(111)
	ax.scatter(curvs,kappas_app-kappas,c=colorset,s=30,lw=0)
	ax.errorbar(curvs_sorted,means,yerr=stddevs,fmt='ko-',lw=1.5,capthick=2)
	ax.grid(True)
	ax.set_xlabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=fsaxlabel)
	ax.set_ylabel(r'$\mathrm{(\mathbf{\kappa}_{app}-\mathbf{\kappa}_{in})\:k_BT}$',fontsize=fsaxlabel)
	ax.set_title('curvature softening',fontsize=fsaxtitle)
	xlims,ylims = ax.get_xlim(),ax.get_ylim()
	ax.set_aspect(abs((xlims[1]-xlims[0])/(ylims[1]-ylims[0]))/1.0)
	plt.savefig(pickles+'fig-meta-meso-kappa-curvature-summary.png',dpi=300)
	if plotshow: plt.show()
	fig.clf()
	



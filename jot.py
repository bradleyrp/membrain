#!/usr/bin/python

thresh = 1

if 0:
	spec_query = [0]
	for s in range(len(spec_query)):
		#---original residuals
		a0 = msets[1].lenscale
		dat0 =array([i for i in mscs[0].t1d[spec_query[s]] if i[0] != 0 and i[0] < thresh])
		dat1 = array([i for i in mscs[1].t1d[spec_query[s]] if i[0] != 0 and i[0] < thresh])
		dat0log = log10(dat0)
		dat1log = log10(dat1)
		cd = scipy.spatial.distance.cdist(dat0log,dat1log)
		cd[cd.argmin(axis=argmax(shape(cd)))]
		resid = mean([cd[i,argmax(shape(cd))] for i in cd.argmin(axis=argmax(shape(cd)))])
		#---new residual signature curves
		datlogs = []
		for m in range(len(mscs)):
			xdat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
			ydat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
			datlogs.append([xdat0,ydat0])
		#---save and report
		print 'result: C_0 = '+('{0:.3f}'.format(c0ask*msets[1].lenscale)).rjust(5)+' (nm^-1)   '+\
			'resid = '+('{0:.3f}'.format(resid)).rjust(10)+'  status: npts = ('+str(len(dat0))+\
			','+str(len(dat1))+')'
		if 'collected_residuals' in globals(): 
			collected_residuals[cgmd_avail.index(cgmd_reference)][s].append([c0ask,resid])
			collected_residual_sigs[cgmd_avail.index(cgmd_reference)][s].append([c0ask,datlogs])


if 0:
	#---codeblock which compares undulations to their residuals, requires mscs
	ax1 = plt.subplot(211)
	ax1.plot(dat0[:,0],dat0[:,1],'r+')
	ax1.plot(dat1[:,0],dat1[:,1],'bx')
	m = 1
	xdat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
	ydat0 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
	ax.plot(xdat0,ydat0,'ro')
	m = 0
	xdat1 = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
	ydat1 = collapse_spectrum(mscs[m].qmagst,mscs[m].t2d[0])
	ax1.plot(xdat1,ydat1,'bo')
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax2 = plt.subplot(212)
	dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
	resid_distn = cd.min(axis=argmax(shape(cd)))
	ax2.plot(dat0log[:,0],dat0log[:,1]-dat1log[:,1],'o-')
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xscale('log')
	plt.show()
	
if 0:	
	#---codeblock which looks at error signatures assuming same wavevectors
	ax2 = plt.subplot(111)
	dat0log = array([i for i in array([xdat0,log10(ydat0)]).T if i[0] < thresh])
	dat1log = array([i for i in array([xdat1,log10(ydat1)]).T if i[0] < thresh])
	ax2.plot(dat0log[:,0],dat0log[:,1]-dat1log[:,1],'o-')
	ax2.set_xscale('log')
	plt.show()
	
#---analyze error signatures
spec_query = [0]
for s in range(len(spec_query)):
	for m in range(2):
		
	

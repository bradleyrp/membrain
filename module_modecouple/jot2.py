#!/usr/bin/python
if 'hqs_hqps' not in globals():
	print 'loading pickles'
	hqs_hqps = unpickle(pickles+'pkl.tmp0.coupling.pkl')
	hqs_h0qps = unpickle(pickles+'pkl.tmp1.coupling.pkl')
	h0qs_hqps = unpickle(pickles+'pkl.tmp2.coupling.pkl')
	h0qs_h0qps = unpickle(pickles+'pkl.tmp3.coupling.pkl')
	qmags = mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi for j in range(0,n2) 
		for i in range(0,m2)])
	print 'constructing kappa'
	kqqp = [[kqs[(int(u/m2)+int(v/n2))%m2,(u%m2+v%n2)%n2] for v in range(m2*n2)] for u in range(m2*n2)]
	print 'constructing matrix'
	bigmat = ((qmags*qmags)*(qmags*qmags)*hqs_hqps+(qmags*qmags)*hqs_h0qps+\
		(qmags*qmags)*h0qs_hqps+h0qs_h0qps)*kqqp
	for lim in [100,200,500,1000,2000,None]:
		st = time.time()
		fullans = linalg.eig(abs(bigmat)[:lim,:lim])
		if lim != None: del fullans
		print 1./60*(time.time()-st)
if 1:
	if 0: plt.imshow(abs(fullans[1]).T,interpolation='nearest',origin='lower',
		norm=(mpl.colors.LogNorm() if 0 else None),cmap=mpl.cm.binary);plt.show()
	if 0: plt.plot(range(4225),fullans[0]);plt.show();
	#fullans = ms.fullans
	#n2,m2 = 65,65
	#import matplotlib as mpl
	rowcombo = mean([array([[fullans[1][fr][j+(i-1)*m2] for j in range(n2)] for i in range(m2)])*fullans[0][fr]**2 for fr in range(m2*n2)],axis=0)
	plt.imshow(abs(rowcombo).T,interpolation='nearest',origin='lower',norm=(mpl.colors.LogNorm() if 0 else None),cmap=mpl.cm.binary);plt.show()

#!/usr/bin/python

#---LOOP SLOW
#-------------------------------------------------------------------------------------------------------------

def ind2to1(i,j):
	return i*n2+j
def ind1to2(k):
	return (int(k/n2),k%n2)
	
tlen = len(msc.hqs)

def term(i,j,k,l):
	if j == 0 and k == 0 and l == 0: print i
	ans = (dot(qmags[i,j],qmags[i,j])*dot(qmags[k,l],qmags[k,l])*\
		abs(mean([msc.hqs[t][i,j]*msc.hqs[t][k,l] for t in range(len(msc.hqs))])) +\
		-1*dot(qmags[i,j],qmags[i,j])*\
		abs(mean([msc.hqs[t][i,j]*msc.cqs[t][k,l] for t in range(len(msc.hqs))])) +\
		-1*dot(qmags[k,l],qmags[k,l])*\
		abs(mean([msc.cqs[t][i,j]*msc.hqs[t][k,l] for t in range(len(msc.hqs))])) +\
		abs(mean([msc.cqs[t][i,j]*msc.cqs[t][k,l] for t in range(len(msc.cqs))]))) *\
		abs(kqs[(i+k)%m2,(j+l)%n2])
	return ans
	
if 0:
	bigmat = [term(ind1to2(u)[0],ind1to2(u)[1],ind1to2(v)[0],ind1to2(v)[1]) 
		for u in range(m2*n2) for v in range(m2*n2)]
	
def term2(u,v):
	ans = (dot(qmags[int(u/n2),u%n2],qmags[int(u/n2),u%n2])*dot(qmags[int(v/n2),v%n2],qmags[int(v/n2),v%n2])*\
		abs(mean([msc.hqs[t][int(u/n2),u%n2]*msc.hqs[t][int(v/n2),v%n2] for t in range(tlen)])) +\
		-1*dot(qmags[int(u/n2),u%n2],qmags[int(u/n2),u%n2])*\
		abs(mean([msc.hqs[t][int(u/n2),u%n2]*msc.cqs[t][int(v/n2),v%n2] for t in range(tlen)])) +\
		-1*dot(qmags[int(v/n2),v%n2],qmags[int(v/n2),v%n2])*\
		abs(mean([msc.cqs[t][int(u/n2),u%n2]*msc.hqs[t][int(v/n2),v%n2] for t in range(tlen)])) +\
		abs(mean([msc.cqs[t][int(u/n2),u%n2]*msc.cqs[t][int(v/n2),v%n2] for t in range(tlen)]))) *\
		abs(kqs[(int(u/n2)+int(v/n2))%m2,(u%n2+v%n2)%n2])
	return ans
	
if 0: 
	bigmat = [[term2(u,v) for v in range(10)] for u in range(10)]
	
#---FAST VECTORIZED
#-------------------------------------------------------------------------------------------------------------

if 0:
	u,v,t = 1000,1000,100
	val = (dot(qmags[int(u/n2),u%n2],qmags[int(u/n2),u%n2])*dot(qmags[int(v/n2),v%n2],qmags[int(v/n2),v%n2])*\
		abs(msc.hqs[t][int(u/n2),u%n2]*msc.hqs[t][int(v/n2),v%n2]) +\
		-1*dot(qmags[int(u/n2),u%n2],qmags[int(u/n2),u%n2])*\
		abs(msc.hqs[t][int(u/n2),u%n2]*msc.cqs[t][int(v/n2),v%n2]) +\
		-1*dot(qmags[int(v/n2),v%n2],qmags[int(v/n2),v%n2])*\
		abs(msc.cqs[t][int(u/n2),u%n2]*msc.hqs[t][int(v/n2),v%n2]) +\
		abs(msc.cqs[t][int(u/n2),u%n2]*msc.cqs[t][int(v/n2),v%n2])) *\
		abs(kqs[(int(u/n2)+int(v/n2))%m2,(u%n2+v%n2)%n2])

if type(msc.hqs) != numpy.ndarray: msc.hqs = array(msc.hqs)
if type(msc.cqs) != numpy.ndarray: msc.cqs = array(msc.cqs)

if 'hqs_hqps' not in globals():
	print 'computing ensemble averaged coupling terms'
	st = time.time()
	hqs_hqps = [[[] for j in range(m2*n2)] for i in range(m2*n2)]
	hqs_h0qps = [[[] for j in range(m2*n2)] for i in range(m2*n2)]
	h0qs_hqps = [[[] for j in range(m2*n2)] for i in range(m2*n2)]
	h0qs_h0qps = [[[] for j in range(m2*n2)] for i in range(m2*n2)]
	for v in range(int(m2*n2)):
		status('v = '+str(v),start=st,i=v,looplen=m2*n2)
		for u in range(int(m2*n2)):
			#---deprecated slow method
			if 0: hqs_hqps[u][v] = [msc.hqs[t][u1,u2]*msc.hqs[t][v1,v2] for t in range(tlen)]
			v1,v2 = int(v/n2),v%n2
			u1,u2 = int(u/n2),u%n2
			hqs_hqps[u][v] = mean(msc.hqs[:,u1,u2]*msc.hqs[:,v1,v2])
			hqs_h0qps[u][v] = mean(msc.hqs[:,u1,u2]*msc.cqs[:,v1,v2])
			h0qs_hqps[u][v] = mean(msc.cqs[:,u1,u2]*msc.hqs[:,v1,v2])
			h0qs_h0qps[u][v] = mean(msc.cqs[:,u1,u2]*msc.cqs[:,v1,v2])
	print 1.*(time.time()-st)/60.
	print 'pickling'
	#pickledump([hqs_hqps,hqs_h0qps,h0qs_hqps,h0qs_h0qps],pickles+'pkl.tmp.coupling.pkl')
	print 'pickled'
	
	if type(hqs_hqps) != numpy.ndarray: hqs_hqps = array(hqs_hqps)
	print 'pickling'
	pickledump(hqs_hqps,pickles+'pkl.v550.tmp0.coupling.pkl')

	if type(h0qs_hqps) != numpy.ndarray: h0qs_hqps = array(h0qs_hqps)
	print 'pickling'
	pickledump(h0qs_hqps,pickles+'pkl.v550.tmp1.coupling.pkl')

	if type(hqs_h0qps) != numpy.ndarray: hqs_h0qps = array(hqs_h0qps)
	print 'pickling'
	pickledump(hqs_h0qps,pickles+'pkl.v550.tmp2.coupling.pkl')

	if type(h0qs_h0qps) != numpy.ndarray: h0qs_h0qps = array(h0qs_h0qps)
	print 'pickling'
	pickledump(h0qs_h0qps,pickles+'pkl.v550.tmp3.coupling.pkl')

	

if 0:
	for lim in [100,200,500,1000,2000,None]:
		st = time.time()
		fullans = linalg.eig(abs(array(hqs_hqps))[:lim,:lim])
		print 1./60*(time.time()-st)
if 0:
	qmags = mset.lenscale*array([ (i-cm)/((Lx)/1.)*2*pi+1j*(j-cn)/((Ly)/1.)*2*pi for j in range(0,n2) for i in range(0,m2)])
	bigmat = (qmags*qmags)*(qmags*qmags)*hqs_hqps+(qmags*qmags)*hqs_h0qps+(qmags*qmags)*h0qs_hqps+h0qs_h0qps
if 0:
	for lim in [100,200,500,1000,2000,None]:
		st = time.time()
		fullans = linalg.eig(abs(bigmat)[:lim,:lim])
		print 1./60*(time.time()-st)
		

#!/usr/bin/python

###############################################
# failed attempt to do fitting to get C0hq
# but do not delete this needs finished
###############################################



from scipy.optimize import leastsq

#---development code
if 0:
	#---generate a simple rule for how C0hq should scale with wavevector
	ref_c0hq_1d = mscs[2].t1d[1]
	ref_c0hq_1d_log = log(ref_c0hq_1d[(ref_c0hq_1d[:,0]>0)])
	[bz,az] = polyfit(ref_c0hq_1d_log[ref_c0hq_1d_log[:,0]<-0.5,0][:len(ref_c0hq_1d_log)/2],
		ref_c0hq_1d_log[ref_c0hq_1d_log[:,0]<-0.5,1][:len(ref_c0hq_1d_log)/2],1)

	#---check the fit
	if 0:
		plt.plot(ref_c0hq_1d[:,0],ref_c0hq_1d[:,1],'bo-')
		plt.plot(ref_c0hq_1d[:,0],[exp(az)*i**(bz) for i in ref_c0hq_1d[:,0]],'ro-')
		plt.xscale('log')
		plt.yscale('log')
		plt.show()

	#---make a dummy hypothesis
	mset = msets[index_md]
	vecs = mean(mset.vecs,axis=0)
	m,n = [i-1 for i in mset.griddims]
	getgrid = array([[[i,j] for j in linspace(0,vecs[1],n)] for i in linspace(0,vecs[0],m)])
	c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
	correlate = (autocorr([c0hypo])*autocorr([msets[0].surf[0][:-1,:-1]]))[0]

	#---convert dummy hypothesis to 1D
	mt,nt = shape(correlate)
	test_c0hq_1d = array([[mscs[index_md].qmagst[i,j],correlate[i,j]] for j in range(nt) for i in range(mt)])

	#---check the plot
	if 0:
		plt.plot(ref_c0hq_1d[:,0],ref_c0hq_1d[:,1],'bo')
		plt.plot(ref_c0hq_1d[:,0],[exp(az)*i**(bz) for i in ref_c0hq_1d[:,0]],'go')
		plt.plot(test_c0hq_1d[:,0],test_c0hq_1d[:,1],'ro')
		plt.xscale('log')
		plt.yscale('log')
		plt.show()
	
	xs = test_c0hq_1d[test_c0hq_1d[:,0]!=0.,0]

	def curvfieldfitter(params,xvals,testsurf ):
		c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
		correlate = (autocorr([c0hypo])*autocorr([msets[0].surf[0][:-1,:-1]]))[0]
		test_c0hq_1d = array([[mscs[index_md].qmagst[i,j],correlate[i,j]] for j in range(nt) for i in range(mt)])
		return array([exp(az)*i**(bz) for i in test_c0hq_1d[test_c0hq_1d[:,0]!=0.,0]])**2-array(test_c0hq_1d[test_c0hq_1d[:,0]!=0.,1])**2

	#p_opt = leastsq(curvfieldfitter,params,args=xs)
	c0best = array([[gauss2d(p_opt[0],getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])

if 0:
	#---define fitting function
	testsurf = msets[0].surf[0]
	def curvfieldfitter(params,xvals):
		c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
		correlate = (autocorr([c0hypo])*autocorr([testsurf[:-1,:-1]]))[0]
		test_c0hq_1d = array([[mscs[index_md].qmagst[i,j],correlate[i,j]] for j in range(nt) for i in range(mt)])
		return array([exp(az)*i**(bz) for i in test_c0hq_1d[test_c0hq_1d[:,0]!=0.,0]])**2-array(test_c0hq_1d[test_c0hq_1d[:,0]!=0.,1])**2

	#---generate a simple rule for how C0hq should scale with wavevector
	ref_c0hq_1d = mscs[2].t1d[1]
	ref_c0hq_1d_log = log(ref_c0hq_1d[(ref_c0hq_1d[:,0]>0)])
	[bz,az] = polyfit(ref_c0hq_1d_log[ref_c0hq_1d_log[:,0]<-0.5,0][:len(ref_c0hq_1d_log)/2],
		ref_c0hq_1d_log[ref_c0hq_1d_log[:,0]<-0.5,1][:len(ref_c0hq_1d_log)/2],1)
		
	popts = []
	for fr in range(len(mset.surf)):
		print fr
		testsurf = mset.surf[fr]
		p_opt = leastsq(curvfieldfitter,params,args=xs)
		print p_opt[0]
		popts.append(p_opt[0])
if 1:
	c0best = array([[gauss2d(mean(popts,axis=0),getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
	plt.imshow(array(c0best).T,interpolation='nearest',origin='lower')
	plt.show()
	besthypo = list(mean(popts,axis=0))[1:2]+[mean(popts,axis=0)[2:4][i]/vecs[i] for i in range(2)]+list(mean(popts,axis=0))[4:]

'''
	
#---compute residual
resid = norm(array([array([exp(az)*i**(bz) for i in test_c0hq_1d[test_c0hq_1d[:,0]!=0.,0]])**2])-[array(test_c0hq_1d[test_c0hq_1d[:,0]!=0.,1])**2],axis=0)


def curvfieldfitter(params):
	c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
	correlate = (autocorr([c0hypo])*autocorr([msets[0].surf[0][:-1,:-1]]))[0]
	test_c0hq_1d = array([[mscs[index_md].qmagst[i,j],correlate[i,j]] for j in range(nt) for i in range(mt)])
	return sum(norm(array([array([exp(az)*i**(bz) 
		for i in test_c0hq_1d[test_c0hq_1d[:,0]!=0.,0]])**2])-
		[array(test_c0hq_1d[test_c0hq_1d[:,0]!=0.,1])**2],axis=0))
	

p_opt = leastsq(curvfieldfitter,params)
#---calculate residual
'''
'''
#def tmpfunc(guess_hypo,ref):
init_hypo = [0.005,2./5,2./5,10,10,0]	
tmp = autocorr([collect_c0s[0][0]])*autocorr([msets[0].surf[0]])[0]
#p_opt = leastsq(gauss2d_residual,array(initpts),args=(target[:,0],target[:,1],target[:,2]))
#params = [0,hypo[0]/10.,vecs[0]*hypo[1],vecs[1]*hypo[2],hypo[3]*10.,hypo[4]*10.,hypo[5]]
#c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
'''

'''
steps
get canonical C0hq curve as a reference
make hypothetical C0 field as start points
compute fft of c0*hq
convert the answer into 1D coordinates
'''

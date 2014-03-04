#!/usr/bin/python

ref_c0hq_1d = mscs[2].t1d[1]
plt.scatter(log10(ref_c0hq_1d[ref_c0hq_1d[:,0]<0.5,0]),log10(ref_c0hq_1d[ref_c0hq_1d[:,0]<0.5,1]));plt.show()
[bz,az] = polyfit(log10(ref_c0hq_1d[ref_c0hq_1d[:,0]<0.5,0]),log10(ref_c0hq_1d[ref_c0hq_1d[:,0]<0.5,1]),1)

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

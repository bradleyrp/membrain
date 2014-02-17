#!/usr/bin/python

def gauss2d(params,x,y):
    '''Two-dimensional Gaussian height function with fluid axis of rotation.'''
    z0,c0,x0,y0,sx,sy,th = params
    #---fix the height shift
    z0 = 0
    return z0+c0*exp(-((x-x0)*cos(th)+(y-y0)*sin(th))**2/2./sx**2)*exp(-(-(x-x0)*sin(th)+
        (y-y0)*cos(th))**2/2./sy**2)

vecs = mean(mset.vecs,axis=0)
m,n = msetmd.griddims
getgrid = array([[[i,j] for j in linspace(0,vecs[1],n)] for i in linspace(0,vecs[0],m)])
params = [0,0.05,vecs[0]/2.,vecs[1]/2.,5,5,0]
c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
	for j in range(n)] for i in range(m)])
plt.imshow(c0hypo.T,interpolation='nearest',origin='lower')
plt.show()


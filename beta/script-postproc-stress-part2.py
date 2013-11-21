#!/usr/bin/python

#---Curvature calculation via stress tensor, post-post processing
#-------------------------------------------------------------------------------------------------------------

#---Parameters
nnnum = 5 #---number of nearest neighbors
numgridpts = 20 #---how many grid points, per direction, to sample (20 is optimal)
exag = -0 #---exaggerates the peaks on the plot
distance_factor = 1 #---exponent on distance scaling in pseudo-RBF
#---Pseudo-Radial basis function computation with K-nearest neighbors
pts = mset.wrappbc(mset.unzipgrid(results/10**exag,vecs=mset.vec(0)),
	vecs=mset.vec(0),mode='grow',growsize=1.5)		
newpts = array([[i,j] for i in linspace(0,vecs[0],numgridpts) for j in linspace(0,vecs[1],numgridpts)])
#---Distance matrix
tree = scipy.spatial.cKDTree(pts[:,0:2])
#---Unpack the points and perform the pseudo-RBF
smoothed = []
for pt in newpts:
	tmp = tree.query(pt,nnnum)
	smoothed.append([pt[0],pt[1],mean([pts[tmp[1][i],2]*1/((tmp[0][i])*sum(tmp[0]))**distance_factor 
		for i in range(len(tmp[0]))])])
nearset = [i for i in smoothed[1:] if np.linalg.norm(i[0:2]-mean(mset.protein[0],axis=0)[0:2])<300.]
#---Drop infinite values and the edges, which are noisy
nearset = [i for i in smoothed if np.isinf(i[2]) == False 
	and (i[0] != 0.0 and i[1] != 0.0 and i[0] != vecs[0] and i[1] != vecs[1])]
meshplot(nearset,show='surf')
meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,color=(0,0,0),opacity=0.5)
#---Print maximum and mean curvatures assuming kappa = 20 kBT
print (max(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
print (min(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
print (mean(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)

if 0:
	#---plot average structure
	meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,color=(0,0,0),opacity=0.5)
	meshplot(mset.surf_mean,vecs=mset.vec(0),show='surf')


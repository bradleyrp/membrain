#!/usr/bin/python

from mayavi import mlab

'''
def normalize_v3(arr):
    #Normalize a numpy array of 3 component vectors shape=(n,3)
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr
'''
if 0:
	vecs = mset.vec(0)
	topxyz_wrapped = self.wrappbc(topxyz,vecs,mode='grow')
	dtri = scipy.spatial.Delaunay(topxyz_wrapped[:,0:2])
	find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices if x in simplex for indx in simplex if indx !=x))
if 0:
	#verts = find_neighbors(0,dtri)
	#pts = [topxyz[i] for i in dtri.vertices[verts[0]]]
	#pts2 = array([topxyz_wrapped[i] for i in dtri.vertices[j] for j in find_neighbors(0,dtri)])
	#pts = array([[topxyz_wrapped[i] for i in dtri.vertices[j]] for j in find_neighbors(0,dtri)])
	#thatwaswrong
	#[list(dtri.vertices[j]) for j in find_neighbors(0,dtri)]

	vertices = topxyz_wrapped
	faces = dtri.simplices
	print 1
	#Create a zeroed array with the same type and shape as our vertices i.e., per vertex normal
	norm = numpy.zeros( vertices.shape, dtype=vertices.dtype )
	#Create an indexed view into the vertex array using the array of three indices for triangles
	tris = vertices[faces]
	#Calculate the normal for all the triangles, by taking the cross product of the vectors v1-v0, and v2-v0 in each triangle             
	n = numpy.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )# n is now an array of normals per triangle. The length of each normal is dependent the vertices, # we need to normalize these, so that our next step weights each normal equally.normalize_v3(n)
	# now we have a normalized array of normals, one per triangle, i.e., per triangle normals.
	# But instead of one per triangle (i.e., flat shading), we add to each vertex in that triangle, 
	# the triangles' normal. Multiple triangles would then contribute to every vertex, so we need to normalize again afterwards.
	# The cool part, we can actually add the normals through an indexed view of our (zeroed) per vertex normal array
	norm[ faces[:,0] ] += n
	norm[ faces[:,1] ] += n
	norm[ faces[:,2] ] += n
	ans = normalize_v3(norm)
	#ans = norm
	toplayer = [ans[i]+topxyz_wrapped[i] for i in range(len(ans))]

if 0:
	meshplot(topxyz,show='wire')
	for i in range(1000):
		print i
		mlab.plot3d(array([topxyz_wrapped[i,0],topxyz_wrapped[i,0]+ans[i,0]]),
			array([topxyz_wrapped[i,1],topxyz_wrapped[i,1]+ans[i,1]]),
			array([topxyz_wrapped[i,2],topxyz_wrapped[i,2]+ans[i,2]]),tube_radius=0.5)
	

'''
procedure
	take a point
	find the neighbor points
	find triangle indices
	
'''

#---restart the effor to code vectors

if 0:
	vecs = mset.vec(0)
	topxyz_wrapped = self.wrappbc(topxyz,vecs,mode='grow')
	dtri = scipy.spatial.Delaunay(topxyz_wrapped[:,0:2])
	find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices if x in simplex for indx in simplex if indx != x))
'''
#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)
'''
if 0:
	# unpacked the find_neighbors and then tried it on a single point to make sure it works
	#meshplot(array([topxyz_wrapped[k] for k in [i for j in [simplex for simplex in dtri.simplices if x in simplex] for i in j]]))
	# so the flattened list of neighboring points in a simplex is:
	#[i for j in [simplex for simplex in dtri.simplices if x in simplex] for i in j]
	# or in list form
	#list(set([i for j in [simplex for simplex in dtri.simplices if x in simplex] for i in j]))
	# they appear to match
	# therefore we should be able to plot one neighborhood as follows
	#meshplot(array([topxyz_wrapped[k] for k in find_neighbors(0,dtri)+[0]]))
	# this works. now let's calculate the normal vectors for each triangle
	# this is where find_neighbors is counterproductive. we really want the triangles separately
	# find_nsimps = lambda x,triang: [list(simplex) for simplex in triang.simplices if x in simplex]
	find_nsimps = lambda x,triang: [list([i for i in simplex if i != x]) for simplex in triang.simplices if x in simplex]
	vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
	# construct the cross products
	#tmp = [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)[0:1]]
	#np.cross(tmp[0][0],tmp[0][1])
	# vectorize these
	#[topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]
	#meshplot(array([topxyz_wrapped[k] for k in find_neighbors(0,dtri)+[0]]))
	#meshpoints([topxyz_wrapped[0]+vecnorm(np.cross(j[0],j[1])) for j in [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]])
	#[topxyz_wrapped[0]+vecnorm(np.cross(j[0],j[1])) for j in [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]]
	# some have the wrong directions
	# to fix this we need to invert the z values for the ones that are in the wrong direction
	# this is a darboux frame hack
	#tmp = [vecnorm(np.cross(j[0],j[1])) for j in [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]]	
	#meshplot(array([topxyz_wrapped[k] for k in find_neighbors(0,dtri)+[0]]))
	#meshpoints([topxyz_wrapped[0]+([i[0],i[1],i[2]] if (i[2] > 0.) else [i[0],i[1],-i[2]]) for i in tmp])
	#meshplot(topxyz_wrapped,show='wire')
	# now we do area weighting
	# tmp here is the normal vectors
	#tmp = [vecnorm(np.cross(j[0],j[1])) for j in [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]]
	# then we need the areas of the triangles
	#tmp2 = [poly_area([topxyz_wrapped[j] for j in k]) for k in [[0]+i for i in find_nsimps(0,dtri)]]
	# then we normalize
	#tmp3 = list(tmp2/sum(tmp2))
	tmp4 = [array(tmp[0])*tmp3[0] for i in range(len(tmp))]
	tmp5 = np.sum(tmp4,axis=0)
	meshplot(array([topxyz_wrapped[k] for k in find_neighbors(0,dtri)+[0]]))
	meshpoints(5.*array(vecnorm(tmp5))+topxyz_wrapped[0])
	#meshplot(topxyz_wrapped,show='wire')
'''	
def surface_norms_slooooooooooow(dtri):
	ans = []
	for ind in range(len(topxyz)):
		print ind
		tmp = [vecnorm(np.cross(j[0],j[1])) for j in [topxyz_wrapped[0]-topxyz_wrapped[i] for i in find_nsimps(0,dtri)]]
		tmp2 = [poly_area([topxyz_wrapped[j] for j in k]) for k in [[0]+i for i in find_nsimps(0,dtri)]]
		tmp3 = list(tmp2/sum(tmp2))
		tmp4 = [array(tmp[0])*tmp3[0] for i in range(len(tmp))]
		tmp5 = np.sum(tmp4,axis=0)
		ans.append(5.*array(vecnorm(tmp5))+topxyz_wrapped[ind])
	return ans
'''	
#def surface_norms():
	# pre-calculate triangle areas
	# pre-compute all triangle normals
	# for each vertex, average them
	
if 0:
	vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
	### new plan	
	# pre-calculate triangle areas
	# pre-compute all triangle normals
	# for each vertex, average them

	pts = topxyz_wrapped
	starttime = time.time()
	#print 'calculating simplex areas'
	#simp_areas = [poly_area([pts[j] for j in i]) for i in dtri.simplices]
	
	# rotating through vertices in a triangle
	ttt = [[(i+j)%3 for i in range(3)] for j in range(3)]
	print 'calculating triangle face normals'
	#trifaces = [[np.cross(pts[j[i[0]]],pts[j[i[1]]]) for i in ttt] for j in dtri.simplices]
	trifaces = [[np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]) for i in ttt] for j in dtri.simplices]
	#needs fixed
	print 'calculating simplex areas'
	simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
	print 1./60.*(time.time()-starttime)	
	#tmp = [array(trifaces[i])*simp_areas[i] for i in range(len(simp_areas))]
	ptsareas = np.zeros(len(dtri.points))
	ptsnorms = np.zeros([len(dtri.points),3])
	print 'summing'
	for s in range(len(dtri.simplices)):
		for p in range(3):
			ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
			#ptsareas[dtri.simplices[s][p]] += simp_areas[s]
	ptsnorms = [2*array(vecnorm(i)) for i in ptsnorms]
	print 1./60.*(time.time()-starttime)
	meshpoints(array(ptsnorms)+array(pts))
	meshplot(topxyz_wrapped)	

if 0:
	# test a single point
	ind = 1
	exneisimps = find_neighbors(ind,dtri)+[ind]
	exnei = array([topxyz_wrapped[k] for k in find_neighbors(ind,dtri)+[ind]])
	#meshplot(array([topxyz_wrapped[k] for k in find_neighbors(0,dtri)+[0]]))
	
	whichsimps = [i for i in range(len(dtri.simplices)) if ind in dtri.simplices[i]]
	j = dtri.simplices[whichsimps[0]]
	tmp = [array(vecnorm(np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]))) for i in ttt]
	
	meshplot(array([topxyz_wrapped[k] for k in find_neighbors(ind,dtri)+[ind]]),show='wire')
	#meshpoints([pts[i]+ptsnorms[(-i+2)%3] for i in list(dtri.simplices[whichsimps[0]])])
	#meshpoints([pts[dtri.simplices[whichsimps[0]][i]]+tmp[i] for i in range(3)],color=(0,0,0))
	#meshpoints(pts[dtri.simplices[whichsimps[0]][2]]+tmp[1],color=(1,1,1))
	sometri = [pts[dtri.simplices[whichsimps[0]][i]] for i in range(3)]
	meshplot(sometri)
	meshpoints([array(vecnorm(np.cross(sometri[1]-sometri[0],sometri[2]-sometri[0]))+sometri[0]),
		array(vecnorm(np.cross(sometri[2]-sometri[1],sometri[0]-sometri[1]))+sometri[1]),
		array(vecnorm(np.cross(sometri[0]-sometri[2],sometri[1]-sometri[2]))+sometri[2])],color=(1,1,1))
	reorder = [0,1,2]
	meshpoints([tmp[i]+pts[j[reorder[i]]] for i in range(3)],color=(0,0,0))
	#meshpoints([tmp[i]+sometri[i] for i in range(3)],color=(0,0,0))
	# checking 90 degrees
	'''
	ex1 = array(vecnorm(np.cross(sometri[1]-sometri[0],sometri[2]-sometri[0]))+sometri[0])
	ex1 = np.cross(sometri[1]-sometri[0],sometri[2]-sometri[0])
	ex2 = sometri[1]-sometri[0]
	ex1 = sometri[2]-sometri[1]
	print np.arccos(np.dot(ex1,ex2)/np.linalg.norm(ex1)/np.linalg.norm(ex2))
	meshpoints(sometri[0]+(ex2))
	meshpoints(sometri[0]+(ex1),color=(1,1,1))
	meshpoints(sometri)
	'''

#---single lipid tilt and density frame
if 0:
	starttime = time.time()
	vecnorm = lambda vec: [i/np.linalg.norm(vec) for i in vec]
	find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices 
		if x in simplex for indx in simplex if indx !=x))
	print 'getting points'
	topxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])
	#botxyz = array([mean(self.universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
	#	for i in self.monolayer_residues[1]])
	toptailxyz = array([mean(self.universe.residues[i].selectAtoms(''.join([i+' or ' for i in director[1:-1]]+[director[-1]])).coordinates(),axis=0) 
		for i in self.monolayer_residues[0]])
	vecs = mset.vec(0)
	topxyz_wrapped = self.wrappbc(topxyz,vecs,mode='grow')
	dtri = scipy.spatial.Delaunay(topxyz_wrapped[:,0:2])
	pts = topxyz_wrapped
	point_permute = [[(i+j)%3 for i in range(3)] for j in range(3)]
	print 'calculating triangle face normals'
	trifaces = [[np.cross(pts[j[i[1]]]-pts[j[i[0]]],pts[j[i[2]]]-pts[j[i[0]]]) 
		for i in point_permute] for j in dtri.simplices]
	print 'calculating simplex areas'
	simp_areas = [abs(1./2*np.dot(vecnorm(i[0]),np.sum(i,axis=0))) for i in trifaces]
	print 1./60.*(time.time()-starttime)	
	ptsareas = np.zeros(len(dtri.points))
	ptsnorms = np.zeros([len(dtri.points),3])
	print 'summing'
	for s in range(len(dtri.simplices)):
		for p in range(3):
			ptsnorms[dtri.simplices[s][p]] += array(vecnorm(trifaces[s][p]))*([1,1,1] 
				if vecnorm(trifaces[s][p])[2] > 0. else [-1,-1,-1])*simp_areas[s]
	ptsnorms = [array(vecnorm(i)) for i in ptsnorms]
if 1:
	print 'calculating angle'
	vecslipids = [toptailxyz[i]-topxyz[i] for i in range(len(topxyz))]
	#print 1./60.*(time.time()-starttime)
	#meshpoints(array(ptsnorms)+array(pts))
	#meshplot(topxyz_wrapped)
	angles = [1./pi*arccos(np.dot(vecslipids[i],ptsnorms[i])/np.linalg.norm(vecslipids[i])/np.linalg.norm(ptsnorms[i])) for i in range(len(topxyz))]
	plotthat = [[topxyz[i][0],topxyz[i][1],50*angles[i]] for i in range(len(topxyz))]
	meshplot(plotthat)


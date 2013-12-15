#!/usr/bin/python

from mayavi import mlab

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

if 1:
	vecs = mset.vec(0)
	topxyz_wrapped = self.wrappbc(topxyz,vecs,mode='grow')
	dtri = scipy.spatial.Delaunay(topxyz_wrapped[:,0:2])
	find_neighbors = lambda x,triang: list(set(indx for simplex in triang.simplices if x in simplex for indx in simplex if indx !=x))

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

if 1:
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

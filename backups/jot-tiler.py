#!/usr/bin/python

tilesize=5
for span in [6]:
	mset.calculate_average_surface()
	vecs = mean(mset.vecs,axis=0)
	#meshplot(mset.surf_mean,vecs=vecs)
	griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
	results2 = numpy.zeros(2*span)
	xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],griddims[1]+2)])
	unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
	interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],fill_value=0.0)
	blocksize = (vecs/griddims)
	#avginterpsurf = array([[i[0]/blocksize[0],i[1]/blocksize[0],(interp(i[0],i[1])+mean(mset.surf_position))/blocksize[2]] for i in xypts])
	avginterpsurf = array([[(interp(i,j)+mean(mset.surf_position))/blocksize[2] for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
	#z0shift = mean(mset.surf_mean+[0,0,mean(mset.surf_position)])/((vecs/griddims)[2])
	#span = spannum*((vecs/griddims)[2]+10**-4)
	#centervblock = (vecs/griddims)[2]
	for pt in dat3dpp:
		if (abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and (abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span):
			results2[pt[2]-avginterpsurf[pt[0]][pt[1]]+span] += (1./2*(pt[3]+pt[7])-pt[11]) 
	#meshplot(results/10**6.+100+40*span,vecs=vecs,opacity=0.5)
	#meshplot(results/10**3.,vecs=vecs,opacity=1.,wirecolor=(1,1,1),lsize=1,show='surf')
	#blurred = array([[mean(results[tilesize*i:(i+1)*tilesize,tilesize*j:(j+1)*tilesize]) for j in range(shape(results)[1]/tilesize)] for i in range(shape(results)[0]/tilesize)])
	#meshplot(blurred/10**1.5,vecs=vecs,show='surf')
	print results[10][10]

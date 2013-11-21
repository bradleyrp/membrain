#!/usr/bin/python -i

if 0:
	from membrainrunner import *
	execfile('plotter.py')
	import scipy.interpolate

	location = 'light'
	execfile('locations.py')

	voxelsize=1.0
	tilesize=5
	mset = unpickle(pickles+'pkl.avgstruct.membrane-v614.md.part0002.skip10.trr.pkl')
	#stress_data_location = '/store-delta/compbio/membrane-v614-enthx4-12800/a4-stress-1.0-framewise/a4a-stress-tensor-framewise/results/'
	stress_data_location = basedir+'/pickle-repository/membrane-v614-stress/'

if 0:
	results_stack = []
	for frame in range(len(mset.surf)):
		print 'Running frame '+str(frame)
		file3dpp = stress_data_location+'md.part0002.fr'+str('%04d'%frame)+'.lp.dat3d'
		dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
		for span in [10]:
			mset.calculate_average_surface()
			vecs = mean(mset.vecs,axis=0)
			griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
			results = numpy.zeros((griddims[0]+1,griddims[1]+1))
			xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],griddims[1]+2)])
			unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
			interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],fill_value=0.0)
			blocksize = (vecs/griddims)
			#avginterpsurf = array([[(interp(i,j)+mean(mset.surf_position))/blocksize[2] for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
			avginterpsurf = array([[(interp(i,j)+mean(mset.surf[frame]))/blocksize[2] for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
			for pt in dat3dpp:
				if (abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and (abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span):
					results[pt[0],pt[1]] += (1./2*(pt[3]+pt[7])-pt[11])*(pt[2]-avginterpsurf[pt[0]][pt[1]])*voxelsize
			#blurred = array([[mean(results[tilesize*i:(i+1)*tilesize,tilesize*j:(j+1)*tilesize]) for j in range(shape(results)[1]/tilesize)] for i in range(shape(results)[0]/tilesize)])
			results_stack.append(results)
			#meshplot(blurred/10**1.5,vecs=vecs,show='surf')
			#raw_input('enter to continue')
			#mlab.clf()

if 0:
	tilesize=5
	blurred_stack = []
	for frame in range(len(results_stack)):
		results = results_stack[frame]
		blurred = array([[mean(results[tilesize*i:(i+1)*tilesize,tilesize*j:(j+1)*tilesize]) for j in range(shape(results)[1]/tilesize)] for i in range(shape(results)[0]/tilesize)])
		blurred_stack.append(blurred)
	avg_result_blurred = mean(blurred_stack,axis=0)
	meshplot(avg_result_blurred/10**3,show='surf')
	#meshpoints(mset.protein[0]+[vecs[0],0,-mean(mset.surf_position)-25],scale_factor=10,color=(1,1,1))
	
if 1:
	tmp=mset.wrappbc(mset.unzipgrid(avg_result/10**3,vecs=mset.vec(0)),vecs=mset.vec(0),mode='grow')
	t0 = time.time()
	#rbfi = scipy.interpolate.Rbf(tmp[:,0],tmp[:,1],tmp[:,2])
	RectBivariateSpline
	interp_obj = scipy.interpolate.RectBivariateSpline(tmp[:,0],tmp[:,1],tmp[:,2])
	print (time.time-t0)/60.
	#xypts = array([[i,j] for i in linspace(0,vecs[0],grid[0]) for j in linspace(0,vecs[1],grid[1])])
	#di = rbfi(xypts[:,0],xypts[:,1])
	#print (time.time-t0)/60.

'''
grid it and plot it
>>> tmp=mset.wrappbc(mset.unzipgrid(avg_result/10**3,vecs=mset.vec(0)),vecs=mset.vec(0),mode='nine')
>>> meshplot(tmp,show='wire')
>>> meshplot(avg_result/10**3,vecs=mset.vec(0))
'''

'''
#OLD PLOTS
if 1: #4-plot
	meshpoints(mset.protein[0]-[0,0,mean(mset.surf_position)+25],scale_factor=10,color=(1,1,1))
	meshpoints(mset.protein[0]+[vecs[0],0,-mean(mset.surf_position)-25],scale_factor=10,color=(1,1,1))
	meshpoints(mset.protein[0]+[-vecs[0],0,-mean(mset.surf_position)-25],scale_factor=10,color=(1,1,1))
	meshpoints(mset.protein[0]+[-2*vecs[0],0,-mean(mset.surf_position)-25],scale_factor=10,color=(1,1,1))
	#meshplot(mset.surf_mean,vecs=vecs)
	checkmesh(mset.surf_mean,vecs=vecs,tess=-1,show='surf')
	meshpoints(mset.protein[0]+[-vecs[0],0,-mean(mset.surf_position)-25],scale_factor=10,color=(1,1,1))
	
if 0: #1 plot
	meshpoints(mset.protein[0]-[0,0,mean(mset.surf_position)+25],scale_factor=10,color=(1,1,1))
	checkmesh(mset.surf_mean,vecs=vecs,tess=-1,show='surf')
'''

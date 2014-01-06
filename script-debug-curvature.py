#!/usr/bin/python -i 

if 0:
	tmp = [result_data.get(['type','target_zones'])[fr][result_data.get(['type','maxhxys'])[fr]] for fr in range(500)]
	meshpoints(mset.protein[0],scale_factor=[20. for i in range(len(mset.protein[0]))])
	tmp2 = [[tmp[i][0],tmp[i][1],10**4*maxhs[i]] for i in range(500)]

	meshpoints(tmp,scale_factor=[20. for i in range(500)],opacity=0.4)

	mset.calculate_average_surface()


	meshplot(mset.surf_mean-30.,vecs=mset.vecs[0],show='surf')
	meshpoints(tmp,scale_factor=[20. for i in range(500)],opacity=1.)

if 0:
	tmp = [result_data.get(['type','target_zones'])[fr][result_data.get(['type','maxhxys'])[fr]] for fr in range(500)]
	tmp2 = [[tmp[i][0],tmp[i][1],maxhs[i]] for i in range(500) if (abs(maxhs[i]) > 0.001 and abs(maxhs[i]) < 0.1)]
	tmp2 = [[tmp[i][0],tmp[i][1],maxhs[i]] for i in range(500)]
	meshpoints(array([i for i in tmp2 if i[2] > 0.]),scale_factor=[abs(10**-1*log(abs(maxhs[i]))) for i in range(len(tmp2)) if tmp2[i][2] > 0.],opacity=1.,color=(1,1,1))
	meshpoints(array([i for i in tmp2 if i[2] < 0.]),scale_factor=[abs(10**-1*log(abs(maxhs[i]))) for i in range(len(tmp2)) if tmp2[i][2] < 0.],opacity=1.,color=(0,0,0))
	meshplot(mset.surf_mean-30.,vecs=mset.vecs[0],show='surf')
	meshpoints(mset.protein[0],scale_factor=[10. for i in range(len(mset.protein[0]))],color=(1,0.8,0.8))

if 1:
	targetxy = [result_data.get(['type','target_zones'])[fr][result_data.get(['type','maxhxys'])[fr]] for fr in range(500)]
	validhispoz = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > curvature_filter[0] and abs(10*maxhs[i]) < curvature_filter[1] and maxhs[i] > 0)]
	validhisneg = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > curvature_filter[0] and abs(10*maxhs[i]) < curvature_filter[1] and maxhs[i] < 0)]
	meshpoints(array([[targetxy[i][0],targetxy[i][1],0.] for i in validhispoz]),scale_factor=[abs(5*10**0*log(abs(maxhs[i]))) for i in validhispoz],opacity=1.,color=(1,1,1))
	meshpoints(array([[targetxy[i][0],targetxy[i][1],0.] for i in validhisneg]),scale_factor=[abs(5*10**0*log(abs(maxhs[i]))) for i in validhisneg],opacity=1.,color=(0,0,0))
	meshplot(mset.surf_mean-30.,vecs=mset.vecs[0],show='surf')
	meshpoints(mset.protein[0],scale_factor=[10. for i in range(len(mset.protein[0]))],color=(1,0.8,0.8))

if 0:
	validhis = [i for i in range(len(maxhs)) if (10*abs(maxhs[i]) > curvature_filter[0] and abs(10*maxhs[i]) < curvature_filter[1])]
	validhs = [10*maxhs[i] for i in validhis]
	fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((1,1),(0,0))
	ax.hist([i for i in validhs],histtype='stepfilled',color=colorcodes[0],alpha=0.8,bins=40,linewidth=2)
	plt.show()

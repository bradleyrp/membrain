#!/usr/bin/python

if 0:
	fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
	target_zones = test.get(['type','target_zones'])
	hmaxraw = []
	center_near_nborhood = [False for i in range(len(test.data))]
	for i in range(len(test.data)):
		if fitted_inds[i]:
			z0,c0,x0,y0,sx,sy,th = test.data[i][0]
			hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
			center_near_nborhood[i] = scipy.spatial.distance.cdist(target_zones[i][:,:2],[[x0,y0]]).min() < 10.
		else: hmaxraw.append(0)
	center_near_nborhood = array(center_near_nborhood)
	magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
		for i in range(len(hmaxraw))])
	cfilt_inds = where(array(1*magfilter+1*center_near_nborhood+1*fitted_inds)==3)[0]
	hmaxdat = array(hmaxraw)[cfilt_inds]
	params = test.get(['type','params'])[cfilt_inds]
	target_zones = test.get(['type','target_zones'])[cfilt_inds]
	resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
		for i in target_zones[j]])) for j in range(len(target_zones))]


if 0:
	vecs = mean(mset.vecs,axis=0)
	fitted_inds = array([type(test.data[i][1]) != list for i in range(len(test.data))])
	target_zones = test.get(['type','target_zones'])
	hmaxraw = []
	center_near_nborhood = [False for i in range(len(test.data))]
	for i in range(len(test.data)):
		if fitted_inds[i]:
			z0,c0,x0,y0,sx,sy,th = test.data[i][0]
			hmaxraw.append(-1*10*curvfac*gauss2dh(test.data[i][0],x0,y0))
			center_near_nborhood[i] = scipy.spatial.distance.cdist(target_zones[i][:,:2],[[x0,y0]]).min() < sqrt(2)*10.
		else: hmaxraw.append(0)
	center_near_nborhood = array(center_near_nborhood)
	magfilter = array([(abs(hmaxraw[i])>smallfilt and abs(hmaxraw[i])<hifilt) \
		for i in range(len(hmaxraw))])
	cfilt_inds = where(array(1*magfilter+1*center_near_nborhood+1*fitted_inds)==3)[0]
	hmaxdat = array(hmaxraw)[cfilt_inds]
	params = test.get(['type','params'])[cfilt_inds]
	target_zones = test.get(['type','target_zones'])[cfilt_inds]
	resids = [sqrt(mean([abs(gauss2d(params[j],i[0],i[1])-i[2])**2 
		for i in target_zones[j]])) for j in range(len(target_zones))]



if 1:
	params = test.get(['type','params'])
	target_zones = test.get(['type','target_zones'])
	thelist = list(where(center_near_nborhood==False)[0])
	#thelist = [thelist[i] for i in [1,4]]
	for i in thelist:
		x0,y0 = params[i][2:4]
		dat = scipy.spatial.distance.cdist(target_zones[i][:,:2],[[x0,y0]])
		j = argmin(dat[:,0])
		#j = list(dat[:,0]).index(sort(dat[:,0])[0])
		meshpoints(array(list(params[i][2:4])+[0.0]),scale_factor=20)
		meshpoints(target_zones[i],scale_factor=10,color=(1,1,1))
		meshpoints(target_zones[i][j],scale_factor=10,color=(1,0,1))
		raw_input('...')
	


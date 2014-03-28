#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

director_aamd_symmetric = ['name P and not resname CHL1','name C218','name C318']
director_aamd_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
	'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = 'name P'

analysis_descriptors = {
	'v509-40000-90000-10':
		{'sysname':'membrane-v509',
		'sysname_lookup':'membrane-v509-ions',
		'trajsel':'s6-kraken-md.part0018.40000-90000-10.ions.xtc',
		'sysname_lookup_struct':'membrane-v509-atomP',
		'trajsel_struct':'s6-kraken-md.part0018.40000-90000-10.atomP.xtc',
		'ionname':'NA',
		'director':director_aamd_symmetric}}
analysis_names = [
	'v509-40000-90000-10'
	][-1:]
routine = ['load','compute','average_leaflets','fit'][3:]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice

if 'load' in routine:
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		trajsel = trajsel_struct
		sysname_lookup = sysname_lookup_struct
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals(),
			keysysname='sysname_lookup_struct',keytrajsel='trajsel_struct')
		mset_surf = MembraneSet()
		mset_surf.load_trajectory((basedir+'/'+grofile,basedir+'/'+trajfile[0]),resolution='aamd')
		#---collect ion positions
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		#---collect monolaye positions
		#---Nb this is where we can refine the method
		mset_surf.identify_monolayers(director)
		surf_select = mset_surf.universe.selectAtoms(selector)
		monoz = []
		for fr in whichframes:
			mset_surf.gotoframe(fr)
			surfpts = surf_select.coordinates()
			monoz.append([mean(surfpts[mset_surf.monolayer_residues[m],2],axis=0) for m in range(2)])

if 'compute' in routine:
      for aname in analysis_names:
		desired_binsize = 1
		upper_binws = [int(round((mset_surf.vec(i)[2]-monoz[i][0])/desired_binsize)) for i in whichframes]
		mid_binws = [int(round((monoz[i][0]-monoz[i][1])/desired_binsize)) for i in whichframes]
		lower_binws = [int(round((monoz[i][1])/desired_binsize)) for i in whichframes]
		#---pick a standard number of bins for consistent zones
		bin_nums = [int(round(mean(i))) for i in [lower_binws,mid_binws,upper_binws]]
		badframes = list(where([ionspos[j][:,2].max()>mset_surf.vec(j)[2] for j in range(5000)])[0])	
		whichframes2 = [i for i in whichframes if i not in badframes]
		binedges = array([np.concatenate((
			linspace(0,monoz[i][1],bin_nums[0])[:-1],
			linspace(monoz[i][1],monoz[i][0],bin_nums[1])[:-1],
			linspace(monoz[i][0],mset_surf.vec(i)[2],bin_nums[2]))) 
			for i in whichframes])
		disctraj  = array([[
			where(ionspos[j][:,2][i]<binedges[j][1:])[0][0] 
			for i in range(len(ionspos[0]))] for j in whichframes2]).T
		bincount = len(binedges[0])-1
		nbins = 100
		jumptimes = [where((disctraj[i,1:]-disctraj[i,:-1])!=0)[0] for i in range(len(disctraj))]
		jumpdurations = [jumptimes[i][1:]-jumptimes[i][:-1] for i in range(len(jumptimes))]
		maxresid = max([max(i) for i in jumpdurations])
		jumphists = [histogram(jumpdurations[i],range=(0,maxresid),bins=nbins)[0] for i in range(len(jumpdurations))]
		jumpbins = [disctraj[i][jumptimes[i][:-1]] for i in range(len(disctraj))]
		residences = [[] for i in range(bincount)]
		for k in range(len(jumpbins)):
			for i in range(len(jumpbins[k])):
				residences[jumpbins[k][i]].append(jumpdurations[k][i])
		maxresid = max([max(i) for i in residences if i != []])
		resdists = [histogram(residences[i],range=(0,maxresid),bins=nbins)[0] for i in range(len(residences)) if i != []]	
		#---accounting
		binlabels = [int(round(i)) for i in (mean(binedges,axis=0)[1:]+mean(binedges,axis=0)[:-1])/2.]
		maxz = binedges[0][-1]
		meanz = mean(mean(monoz,axis=0))
		howfar = [i-meanz-maxz*((i-meanz)>maxz/2.) for i in binlabels]
		'''
		#---plot some single-frame distributions
		for i in range(0,shape(disctraj)[1],1):
		#for i in [1]:
			hist,edges = histogram(disctraj[:,i],range=(0,bincount),bins=bincount);
			mids = array(edges[1:]+edges[:-1])/2.
			plt.plot(mids,hist);
		plt.show()
		'''

		hists = []
		meanbinedges = mean(binedges,axis=0)	
		monobins = [array([abs(mean(monoz,axis=0)[i]-meanbinedges[j])
			for j in range(len(meanbinedges))]).argmin() for i in range(2)][::-1]
		for i in range(0,shape(disctraj)[1],1):
			hist,edges = histogram(disctraj[:,i],range=(0,bincount),bins=bincount);
			hists.append(hist)
		mids = array(edges[1:]+edges[:-1])/2.
		bw = mean(meanbinedges[1:]-meanbinedges[:-1]) # bin width in Angstroms.
		# plt.axvline(x=bw*edges[monobins[0]],ymin=0.,ymax=1.)
		# plt.axvline(x=bw*edges[monobins[1]],ymin=0.,ymax=1.)
		# plt.plot(bw*mids,mean(hists,axis=0),'o-');
		# plt.show()

		print 'Splitting the system into two leaflets'
		zeroes = where(mean(hists,axis=0)==0)
		middle_bin = int(mean(zeroes))
		upper = [i > 0 for i in howfar]
		upper_positions = [i for i in binlabels*array(upper) if i != 0]
		lower_positions = [i for i in binlabels*~array(upper) if i != 0]

		upper_bins = [binlabels.index(upper_positions[j]) for j in range(len(upper_positions))]
		lower_bins = [binlabels.index(lower_positions[j]) for j in range(len(lower_positions))]
		print 'Upper bins: '+str(len(upper_positions))
		print 'Lower bins: '+str(len(lower_positions))
		bincount = len(upper_bins)+len(lower_bins)
		hists = []
		''' This bollocks is hopelessly broken.
		for i in range(0,shape(disctraj)[1],1):
			hist,edges = histogram([j for j in disctraj[:,i] if j in upper_bins],bins=len(upper_bins));
			hists.append(hist)
			plt.plot(upper_bins,hist, 'bo-')
		#plt.plot([j*bw for j in mids[upper_bins]],mean(hists,axis=0),'bo-');
		hists = []        
		for i in range(0,shape(disctraj)[1],1):
			hist,edges = histogram([j for j in disctraj[:,i] if j in lower_bins],bins=len(lower_bins));
			hists.append(hist)
			plt.plot(lower_bins,hist, 'ro-')
		#plt.plot([j*bw for j in mids[lower_bins]],mean(hists,axis=0),'ro-');
		'''
		for i in range(0,shape(disctraj)[1],1):
			hist,edges = histogram(disctraj[:,i],bins=bincount);
			hists.append(hist)
		#plt.plot(mids,mean(hists,axis=0))		
		#plt.plot([bw*mids[i] for i in lower_bins],[mean(hists,axis=0)[i] for i in lower_bins],'ro-')
		#plt.plot([bw*mids[i] for i in upper_bins],[mean(hists,axis=0)[i] for i in upper_bins],'bo-')
		#plt.show()
		
		lower_hist = [mean(hists,axis=0)[i] for i in lower_bins] # Stop computing this a thousand times.
		upper_hist = [mean(hists,axis=0)[i] for i in upper_bins]
		lower_peak_max = argmax(lower_hist)
		lower_peak_pos = bw*mids[lower_bins[lower_peak_max]]
		upper_peak_max = argmax(upper_hist)
		upper_peak_pos = bw*mids[upper_bins[upper_peak_max]]
		'''
		diff = [lower_bins[i+1] - lower_bins[i] for i in range(len(lower_bins)-1)]
		lower_discontinuity = np.where(array(diff) > 1)
		'''
		tmp = []
		for i in lower_bins:
			old_position = bw*mids[i]
			# Reflect the lower peak
			if i == lower_peak_max:
				tmp.append(lower_peak_pos)
			elif i > upper_bins[upper_peak_max]: # This needs to be wrapped.
				# Move this to start up after what happens when old position is zero.
				tmp.append(2*lower_peak_pos + abs(bw*mids[-1] - old_position))
			elif i > lower_peak_max: # This is effectively inside the bilayer.
				tmp.append(lower_peak_pos - abs(lower_peak_pos - old_position))
			elif i < lower_peak_max: # This is above the bilayer.
				tmp.append(lower_peak_pos + abs(lower_peak_pos - old_position))
			else:
				print "Not enough coffee today to move the peaks."
		# Offset the lower peak
		lower_position = [] # This will contain the final adjusted positions for the "bottom" leaflet
		for i in range(len(tmp)):
			lower_position.append(tmp[i]-lower_peak_pos)
			
		upper_position = [] # This will contain the final adjusted positions for the "top" leaflet
		for i in upper_bins:
			old_position = bw*mids[i]
			# Offset only
			upper_position.append(old_position-upper_peak_pos)
		
		
		#plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet')
		#plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet')
		#plt.show()
		
		lower = array(zip(lower_position, lower_hist))
		upper = array(zip(upper_position, upper_hist))
		lower_sorted = lower[np.argsort(lower[:, 0])]
		upper_sorted = upper[np.argsort(upper[:, 0])]

		# data[:,n] -- get entire column of index n
		# argsort() -- get the indices that would sort it
		# data[data[:,n].argsort()] -- get data array sorted by n-th column
		
		symmetric_hist = []
		symmetric_pos = []
		# This is unpythonic but it works.
		if 'average_leaflets' in routine:
			if len(lower_bins) == len(upper_bins):
				for i in range(len(lower_bins)):
					symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
					symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
			elif len(lower_bins) < len(upper_bins):
				print 'Different number of bins, so double check the averaging!'
				for i in range(len(lower_bins)):
					symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
					symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
				for i in range((len(upper_bins) - len(lower_bins))):
					last_step = len(lower_bins)
					symmetric_hist.append(upper_sorted[last_step+i][1])
					symmetric_pos.append(upper_sorted[last_step+i][0])
			else:
				print 'Different number of bins, so double check the averaging!'
				for i in range(len(upper_bins)):
					symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
					symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
				for i in range((len(lower_bins) - len(upper_sorted))):
					last_step = len(upper_bins)
					symmetric_hist.append(lower_sorted[last_step+i][1])
					symmetric_pos.append(lower_sorted[last_step+i][0])
		else:
			print 'Pick a leaflet to fit and set that equal to symmetric_hist and symmetric_pos'
		#plt.figure()
		#plt.scatter(symmetric_pos,symmetric_hist, c='g', s=80, alpha=0.3, label='Average')
		#plt.legend()
		#plt.show()
	
if 'fit' in routine:
	from scipy.optimize import curve_fit
	def func(z, kappa, phi):
		# result is the ratio of ion at point (z) to bulk value.
		# This is positive for positive ions. For negative ions, flip the signs!
		result = ((1 + exp(-kappa*z) * math.tanh(phi/4.)) / (1 - exp(-kappa*z) * math.tanh(phi/4.)))**2
		return result
			
	guess = [5, 2] # Debye screening should be ~0.5 nm for 1 M sodium in water
	
	# Confusingly, it looks like *more* weight means the points count *less*
	#weights = [1/sqrt(array(symmetric_pos)[array(symmetric_pos)>0][i])\
	#	for i in range(len(array(symmetric_pos)[array(symmetric_pos)>0]))]
	weights = [1.0 for i in range(len(array(symmetric_pos)[array(symmetric_pos)>0]))]

	# Fit positive values only.
	fit, fit_covariance = curve_fit(func, array(symmetric_pos)[array(symmetric_pos)>0] , array(symmetric_hist)[array(symmetric_pos)>0], guess, weights)
	sigma = [sqrt(fit_covariance[0,0]), sqrt(fit_covariance[1,1])]
	
	print 'FIT'
	print 'kappa = '+str(fit[0])+' +/- '+str(sigma[0])
	print 'phi(s) = '+str(fit[1])+' +/- '+str(sigma[1])
	
	y = func(array(symmetric_pos)[array(symmetric_pos)>0],fit[0],fit[1])

	plt.figure()
	plt.scatter(symmetric_pos,symmetric_hist, c='g', s=80, alpha=0.3, label='Average')
	plt.plot(array(symmetric_pos)[array(symmetric_pos)>0], y, c='k', lw=4, label='Fit')
	#plt.plot(array(symmetric_pos), func(array(symmetric_pos), fit[0], fit[1]),\
	#array(symmetric_pos), func(array(symmetric_pos), fit[0] + sigma[0], fit[1] - sigma[1]),\
	#array(symmetric_pos), func(array(symmetric_pos), fit[0] - sigma[0], fit[1] + sigma[1])\
	#	)
	plt.legend()
	plt.show()
		

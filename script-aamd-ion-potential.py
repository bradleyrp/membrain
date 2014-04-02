#!/usr/bin/python


if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

font = {'family' : 'sans-serif',
        'size'   : 22}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}',r'\usepackage{amsmath}',
					r'\usepackage{siunitx}',r'\sisetup{detect-all}',
		                        r'\usepackage{helvet}',r'\usepackage{sansmath}',
		                        r'\sansmath', r'\usepackage{upgreek}']
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['ytick.major.pad'] = 8



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
	 'ion_positive':'NA',
	 'ion_negative':'CL',
	 'director':director_aamd_symmetric}}
analysis_names = [
	'v509-40000-90000-10'
	][-1:]
routine = ['load','compute','average_leaflets','fit_together'][-1:]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

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
		ionspos_pos = []
		ionspos_neg = []
		ion_select_pos = mset.universe.selectAtoms('name '+ion_positive)
		ion_select_neg = mset.universe.selectAtoms('name '+ion_negative)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos_pos.append(ion_select_pos.coordinates())
			ionspos_neg.append(ion_select_neg.coordinates())
			clock.append(mset.universe.trajectory[fr].time)

		mset_surf.identify_monolayers(director)
		surf_select = mset_surf.universe.selectAtoms(selector)
		monoz = []
		for fr in whichframes:
			mset_surf.gotoframe(fr)
			surfpts = surf_select.coordinates()
			monoz.append([mean(surfpts[mset_surf.monolayer_residues[m],2],axis=0) for m in range(2)])

if 'compute' in routine:
	for aname in analysis_names:
		ions = [ion_positive,ion_negative]
		for ionname in ions:
			if str(ionname) == str(ion_positive):
				ionspos = ionspos_pos
				print 'Working on positive ion'
			else: 
				ionspos = ionspos_neg
				print 'Working on negative ion'

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
			'''
			# Now, keep track of the actual distances.
			####################################################################################
			d_to_top_leaf = [bw*(i-j) for i, j in zip (mids,[monoz[k][0] for k in range(len(monoz))])] # For all ions.
			d_top = [d_to_top_leaf[i] for i in upper_bins]
			d_bottom = [(i-j) for i, j in zip (tmp,[monoz[k][1] for k in range(len(monoz))])] # Just for lower ions.
			####################################################################################
			'''
			
			# Offset the lower peak
			if ionname != 'CL' or ionname != 'Cl':
				lower_position = [] # This will contain the final adjusted positions for the "bottom" leaflet
				for i in range(len(tmp)):
					lower_position.append(tmp[i]-lower_peak_pos)
					
			upper_position = [] # This will contain the final adjusted positions for the "top" leaflet
			for i in upper_bins:
				old_position = bw*mids[i]
				# Offset only
				upper_position.append(old_position-upper_peak_pos)
			show_each_leaflet = 0	
			if 'show_each_leaflet':
				plt.title('Peaks aligned?')
				plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet')
				plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet')
				plt.show()

				
			# After the shifting, resort the bins by distance from membrane surface
			lower = array(zip(lower_position, lower_hist))
			upper = array(zip(upper_position, upper_hist))
			lower_sorted = lower[np.argsort(lower[:, 0])]
			upper_sorted = upper[np.argsort(upper[:, 0])]

			if ionname == 'CL' or ionname == 'Cl':
				# In this case, we don't want to shift to the max. We want to shift to the first nonzero element.
				# IOW, shift so that the inflection of the curve is positive (and thus, fitted).
				# This is more easily done after the sorting!
				lower_positive = lower_sorted[:,1]>0
				position_of_rise = lower_sorted[lower_sorted[:,1]>0][0][0]
				diff = 0.0 - position_of_rise		
				for i in range(len(lower_sorted)):
					#lower_sorted[i][0] += diff
					lower_sorted[i][1] *= -1.0
					
				upper_positive = upper_sorted[:,1]>0
				position_of_rise = upper_sorted[upper_sorted[:,1]>0][0][0]
				diff = 0.0 - position_of_rise	
				for i in range(len(upper_sorted)):
					#upper_sorted[i][0] += diff
					upper_sorted[i][1] *= -1.0
					
				show_each_leaflet = 0
				if 'show_each_leaflet':
					plt.title('Is the inflection at positive distance?')
					plt.scatter([lower_sorted[i][0] for i in range(len(lower_sorted))],[lower_sorted[i][1] for i in range(len(lower_sorted))],c='r', s=40, alpha=0.5, label='Lower leaflet')
					plt.scatter([upper_sorted[i][0] for i in range(len(upper_sorted))],[upper_sorted[i][1] for i in range(len(upper_sorted))],c='b', s=40, alpha=0.5, label='Upper leaflet')
					plt.show()				

			'''
			bottom = array(zip(d_bottom, lower_hist))
			top = array(zip(d_top, upper_hist))
			bottom_sorted = bottom[np.argsort(bottom[:, 0])]
			top_sorted = top[np.argsort(top[:, 0])]
			plt.title('Actual distances?')
			plt.scatter([bottom_sorted[i][0] for i in range(len(bottom_sorted))], [lower_sorted[i][1] for i in range(len(lower_sorted))], c='r', s=40, alpha=0.5, label='Lower leaflet')
			plt.scatter([top_sorted[i][0] for i in range(len(top_sorted))], [upper_sorted[i][1] for i in range(len(upper_sorted))] ,c='b', s=40, alpha=0.5, label='Upper leaflet')
			plt.show()
			'''
			###################################################################################
			#raw_upper_position = upper_position+upper_peak_pos
			#real_upper_position = raw_upper_position - mean([monoz[k][0] for k in range(len(monoz))]) 
			# Really should be doing this framewise, but this saves a lot of time and is easier and is okay,
			# because the variance is <0.3 A.
			# numpy.var([monoz[k][0] for k in range(len(monoz))])
			#raw_lower_position = lower_position+lower_peak_pos
			#real_lower_position = raw_lower_position - mean([monoz[k][1] for k in range(len(monoz))]) 
			###################################################################################
			peak_bins = np.array(mean(hists,axis=0)).argsort()[-2:][::-1]
			differences = abs(peak_bins - middle_bin)
			peak_to_center = bw*mean(differences)
			P_to_P_distance = mean([monoz[k][0]-monoz[k][1] for k in range(len(monoz))])
			peak_to_P_surface = peak_to_center - P_to_P_distance/2.
		
			symmetric_hist = []
			symmetric_pos = []
			# This is unpythonic but it works.
			if 'average_leaflets'  in routine and ionname != 'CL':
				if len(lower_bins) == len(upper_bins):
					for i in range(len(lower_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
				elif len(lower_bins) < len(upper_bins):
					print 'Different number of bins, so double check the averaging!'
					print "plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet'); plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet'); plt.show()"
					for i in range(len(lower_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
					for i in range((len(upper_bins) - len(lower_bins))):
						last_step = len(lower_bins)
						symmetric_hist.append(upper_sorted[last_step+i][1])
						symmetric_pos.append(upper_sorted[last_step+i][0])
				else:
					print 'Different number of bins, so double check the averaging!'
					print "plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet'); plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet'); plt.show()"
					for i in range(len(upper_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
					for i in range((len(lower_bins) - len(upper_sorted))):
						last_step = len(upper_bins)
						symmetric_hist.append(lower_sorted[last_step+i][1])
						symmetric_pos.append(lower_sorted[last_step+i][0])
						
				pos_pos = array([i for i in symmetric_pos])
				pos_hist = array([i for i in symmetric_hist])
				pos_peak_to_center = peak_to_center
				pos_peak_to_P_surface = peak_to_P_surface
									
			elif 'average_leaflets' in routine and (ionname == 'CL' or ionname == 'Cl'):
				if len(lower_bins) == len(upper_bins):
					for i in range(len(lower_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
				elif len(lower_bins) < len(upper_bins):
					print 'Different number of bins, so double check the averaging!'
					print "plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet'); plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet'); plt.show()"
					for i in range(len(lower_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
					for i in range((len(upper_bins) - len(lower_bins))):
						last_step = len(lower_bins)
						symmetric_hist.append(upper_sorted[last_step+i][1])
						symmetric_pos.append(upper_sorted[last_step+i][0])
				else:
					print 'Different number of bins, so double check the averaging!'
					print "plt.scatter(lower_position,lower_hist,c='r', s=40, alpha=0.5, label='Lower leaflet'); plt.scatter(upper_position,upper_hist,c='b', s=40, alpha=0.5, label='Upper leaflet'); plt.show()"
					for i in range(len(upper_bins)):
						symmetric_hist.append(mean([lower_sorted[i][1], upper_sorted[i][1]]))
						symmetric_pos.append(mean([lower_sorted[i][0], upper_sorted[i][0]]))
					for i in range((len(lower_bins) - len(upper_sorted))):
						last_step = len(upper_bins)
						symmetric_hist.append(lower_sorted[last_step+i][1])
						symmetric_pos.append(lower_sorted[last_step+i][0])
									
				neg_pos = array([i for i in symmetric_pos])
				neg_hist = array([i for i in symmetric_hist])
				neg_peak_to_center = peak_to_center
				neg_peak_to_P_surface = peak_to_P_surface
			
			elif 'average_leaflets' in routine:
				print "I'm not sure if I should adjust the peaks for a positive or negative ion."
			else:
				print 'Pick a leaflet to fit and set that equal to symmetric_hist and symmetric_pos'
			
if 'fit' in routine:
	from scipy.optimize import curve_fit
	for symmetric_pos in range(2): 
		if symmetric_pos == 0:
			symmetric_pos = pos_pos
			symmetric_hist = pos_hist
			positive = 1
		else:
			symmetric_pos = neg_pos
			symmetric_hist = neg_hist
			positive = -1
			
		def func(z, kappa, phi, z0, bulk):
			# result is the ratio of ion at point (z) to bulk value.
			if positive == 1:
				result = bulk*((1 + exp(-kappa*(z-z0)) * math.tanh(phi/4.)) / \
					       (1 - exp(-kappa*(z-z0)) * math.tanh(phi/4.)))**2
			else:
				result = bulk*((1 - exp(-kappa*(z-z0)) * math.tanh(phi/4.)) / \
					       (1 + exp(-kappa*(z-z0)) * math.tanh(phi/4.)))**2
			return result
			
		guess = [5, 2, 1, 5] # Debye screening should be ~0.5 nm for 1 M sodium in water
		# Confusingly, it looks like *more* weight means the points count *less*
		weights = [1.0 for i in range(len(symmetric_pos[symmetric_pos>0]))]
		weights = array(weights)/sum(weights)

		# Fit distance positive values only.
		fit, fit_covariance = curve_fit(func, symmetric_pos[symmetric_pos>0], \
						symmetric_hist[symmetric_pos>0], array(guess), weights)
		sigma = [sqrt(fit_covariance[0,0]), sqrt(fit_covariance[1,1]), \
			 sqrt(fit_covariance[2,2]), sqrt(fit_covariance[3,3])]
	
		print 'kappa = '+str(fit[0])+' +/- '+str(sigma[0])
		print 'phi(s) = '+str(fit[1])+' +/- '+str(sigma[1])
		print 'z offset = '+str(fit[2])+' +/- '+str(sigma[2])
		print '[bulk] = '+str(fit[3])+' +/- '+str(sigma[3])
	
		y = func(symmetric_pos[symmetric_pos>0],fit[0],fit[1],fit[2],fit[3])

		plt.figure()
		plt.scatter(symmetric_pos,symmetric_hist, c='g', s=80, alpha=0.3, label='Average')
		plt.plot(symmetric_pos[symmetric_pos>0], y, c='k', lw=4, label='Fit')
		'''
		if 'plot_with_uncertainy':
			plt.plot(array(symmetric_pos), func(array(symmetric_pos), fit[0], fit[1]),\
				 array(symmetric_pos), func(array(symmetric_pos), fit[0] + sigma[0], fit[1] - sigma[1]),\
				 array(symmetric_pos), func(array(symmetric_pos), fit[0] - sigma[0], fit[1] + sigma[1]))
		'''		 
		plt.legend()
		plt.show()

def pot_func(params,pos,valence):
	kappa, phi, z0, bulkp, bulkn = params
	#z0 = 0
	resids = zeros(len(valences))
	if valence == 1:
		val = bulkp*((1 + exp(-kappa*(pos-z0)) * math.tanh(phi/4.)) /
			(1 - exp(-kappa*(pos-z0)) * math.tanh(phi/4.)))**2
	elif valence == -1:
		val = bulkn*((1 - exp(-kappa*(pos-z0)) * math.tanh(phi/4.)) /
			(1 + exp(-kappa*(pos-z0)) * math.tanh(phi/4.)))**2
	return val

def pot_resids(params,pos,hist,valence):
	kappa, phi, z0, bulkp, bulkn = params
	#z0 = 0
	resids = zeros(len(valences))
	for i in range(len(valence)):
		if valence[i] == 1:
			resids[i] = (hist[i] - bulkp*((1 + exp(-kappa*(pos[i]-z0)) * math.tanh(phi/4.)) /
				(1 - exp(-kappa*(pos[i]-z0)) * math.tanh(phi/4.)))**2)**2
		elif valence[i] == -1:
			resids[i] = (hist[i] - bulkn*((1 - exp(-kappa*(pos[i]-z0)) * math.tanh(phi/4.)) /
				(1 + exp(-kappa*(pos[i]-z0)) * math.tanh(phi/4.)))**2)**2
	return resids

if 'fit_together' in routine:
	from scipy.optimize import curve_fit
	from scipy.optimize import leastsq

	inds_pos = (pos_pos>0)
	inds_neg = (neg_pos>-12)
	pos_cut = np.concatenate((pos_pos[inds_pos],neg_pos[inds_neg]))
	hist_cut = np.concatenate((pos_hist[inds_pos],neg_hist[inds_neg]))
	valences = array([1 for i in range(sum(inds_pos))]+[-1 for i in range(sum(inds_neg))])

	posst = [0.08177729, 2.30000692, 0., 2.92884602, -4.]
	negst = [0.08177729, 2.30000692, -15., 2.92884602, -4.]
	print 'Fitting negative ion'
	p_opt_neg = leastsq(pot_resids,array(negst),
		args=(pos_cut[sum(inds_pos):],hist_cut[sum(inds_pos):],valences[sum(inds_pos):]))
	print p_opt_neg	
	print 'Fitting positive ion'
	p_opt_pos = leastsq(pot_resids,array(posst),
		args=(pos_cut[2:sum(inds_pos)],hist_cut[2:sum(inds_pos)],valences[2:sum(inds_pos)]))
	print p_opt_pos	

	pos_pos_cut = pos_pos[inds_pos]
	hist_pos_cut = pos_hist[inds_pos]	
	pos_neg_cut = neg_pos[inds_neg]
	hist_neg_cut = neg_hist[inds_neg]
		
	fig = plt.figure(figsize=(11,8.5))
	gs = gridspec.GridSpec(1,2,wspace=0.2,hspace=0.05)
	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1])
	
	#ax1.plot(pos_pos_cut-p_opt_pos[0][2],[pot_func(p_opt_pos[0],pos_pos_cut[i],1)/p_opt_pos[0][3] for i in range(len(pos_pos_cut))],'bo-')
	#ax1.plot(pos_pos_cut,array(hist_pos_cut)/p_opt_pos[0][3],'ro-') # Fitted data
	#ax1.bar(pos_pos,pos_hist/p_opt_pos[0][3],width=bw, color='g',alpha=0.2) # Real data
	
	#ax2.plot(pos_neg_cut-p_opt_neg[0][2],[pot_func(p_opt_neg[0],pos_neg_cut[i],-1)/p_opt_neg[0][4] for i in range(len(pos_neg_cut))],'bo-')
	
	#ax2.plot(pos_neg_cut,array(hist_neg_cut)/p_opt_neg[0][4],'ro-') # Fitted data
	#ax2.bar(neg_pos,neg_hist/p_opt_neg[0][4],width=bw, color='g', alpha=0.2) # Real data
	
	#ax.axvline(x=0, ymin=0, ymax=1, color='k')
	ax1.axvspan(-20, 0, facecolor='0.5', alpha=0.2)
	ax2.axvspan(-20, 0, facecolor='0.5', alpha=0.2)
	#ax.annotate('Membrane surface', xy=(0, 10), xytext=(-40, 10), arrowprops=dict(facecolor='black', shrink=0.05, width=0.5))
	#ax1.annotate('Membrane',  xy=(-20, 10))
	
	pos_real_data_shift = pos_peak_to_P_surface - pos_pos[pos_hist.argmax()]
	neg_real_data_shift = neg_peak_to_P_surface - neg_pos[neg_hist.argmin()]

	ax1.bar(pos_pos+pos_real_data_shift,pos_hist/p_opt_pos[0][3],width=bw, color=clrs[2],alpha=0.2) # Synchronized data!
	ax1.plot(pos_pos_cut+pos_real_data_shift,[pot_func(p_opt_pos[0],pos_pos_cut[i],1)/p_opt_pos[0][3] for i in range(len(pos_pos_cut))],'ko-',label='$\kappa$ = %.2f, $\phi$ = %.2f' %(p_opt_pos[0][0],p_opt_pos[0][1])) # Fit
	ax2.bar(neg_pos+neg_real_data_shift,neg_hist/p_opt_neg[0][4],width=bw, color=clrs[2],alpha=0.2) # Synchronized data!
	ax2.plot(pos_neg_cut+neg_real_data_shift,[pot_func(p_opt_neg[0],pos_neg_cut[i],-1)/p_opt_neg[0][4] for i in range(len(pos_neg_cut))],'ko-',label='$\kappa$ = %.2f, $\phi$ = %.2f ' %(p_opt_neg[0][0],p_opt_neg[0][1])) # Fit

	ax1.legend(fontsize=16)
	ax2.legend(fontsize=16)
	#ax1.set_xlabel('Distance from membrane surface (\AA)',fontsize=22)
	ax1.set_ylabel('Concentration relative to bulk', fontsize=22)
	fig.text(0.25, 0.02, "Distance from membrane surface (\AA)",fontsize=22)
	gs.tight_layout(fig, rect=[0.0, 0.03, 1, 1]) # Leave space for the common x-label.
	plt.savefig(pickles+'fig-'+aname+'-ion-potential.png',dpi=300)
	plt.close()
	#ax.set_yscale('log')
	#plt.show()			
	
	


#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

analysis_descriptors = {
'v530-30000-100000-100':
		{'sysname':'membrane-v530',
		'sysname_lookup':'membrane-v530-ions',
		'trajsel':'u5-sim-trestles-md.part0006.30000-100000-100.ions.xtc',
	'structure_pkl':
			'pkl.structures.membrane-v530.a4-surfacer.u5-sim-trestles-md.part0006.30000-100000-100.pkl',
		'ionname':'NA'}}
analysis_names = ['v530-30000-100000-100']


routine = ['compute','postproc','computexyz',][0:2]

#---method
upto = 500 #---how far to only look at the diffusion curves
timelimit = 3*10**2 #---maximum ps for fitting which depends on occupancy because displacements shrink
occupancy = 1. #---occupancy rate for consideration of the displacement
bwid = 5 #---angstrom width of slices where default is 10 and zonesub must be none
zonesub = [1,2,3,4,5,6,7,8,14,15,16,17,18,19,20,21] #---custom selection of zones

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice
if 'compute' in routine or 'computexyz' in routine or 'mastermsd_zones' not in globals():
	for aname in analysis_names:
		for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
		#---load
		print 'status: loading trajectory'
		grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
		mset_surf = unpickle(pickles+structure_pkl)
		#---no looping over trajfile names, so only specify one in the analysis_descriptors
		traj = trajfile[0]
		mset.load_trajectory((basedir+'/'+grofile,basedir+'/'+traj),resolution='aamd')
		checktime()
		#---check for pre-existing pickle
		resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
		#---disabled file-checking for now
		if 0:
			mdionskate = unpickle(pickles+resultpkl)
			if mdionskate != None:
				raise Exception('except: pkl already exists so figure out your naming problems')
		#---get ion positions and times
		print 'status: getting ions'
		clock = []
		ionspos = []
		ion_select = mset.universe.selectAtoms('name '+ionname)
		whichframes = range(len(mset.universe.trajectory))
		for fr in whichframes:
			mset.gotoframe(fr)
			ionspos.append(ion_select.coordinates())
			clock.append(mset.universe.trajectory[fr].time)
		vecs=mean(mset_surf.vecs,axis=0)
#		ionspos = array(ionspos)[:-1] # Maybe not necessary if the times match?
		print len(ionspos)
		ionstraj = []
		for ind in range(shape(ionspos)[1]):
			course = array(ionspos)[:,ind]
			#---three-line handling PBCs
			hoplistp = (course[1:]-course[:-1])>array(mset_surf.vecs)[1:]/2.
			hoplistn = (course[:-1]-course[1:])>array(mset_surf.vecs)[1:]/2.
			course_nojump = course[1:]-(cumsum(1*hoplistp-1*hoplistn,axis=0))*array(mset_surf.vecs)[1:]
			ionstraj.append(course_nojump)
		ionstraj=array(ionstraj)
		nions = len(ionstraj)
		nframes = len(ionstraj[0])
		#---specify zones
		center = mean(mset_surf.surf_position)
		cslice = int(round(center/bwid))
		thick = mean(mset_surf.surf_thick)
		zonesabs = list(arange(center,0-bwid,-bwid))[::-1][:-1]+list(arange(center,vecs[2]+bwid,bwid))
		zonecenters = array(zonesabs[:-1]+zonesabs[1:])/2.
		rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+\
			(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)]
		roundz = [[int(round(i)) for i in (rewrapz[j]-center)/bwid] for j in range(nions)]
		roundz = roundz - array(roundz).min()
		nzones = ptp(roundz)
		#---updated for subsetting the zones to save memory when using small bin widths
		nzones = len(zonesub)
		#---select subset of ions if desired
		ionsel = slice(0,nions)
		nions = len(range(nions)[ionsel])
		#---if you ask for xyz diffusion, it skips the z-decomposition entirely
		if 'computexyz' in routine:
			print 'status: precomputing displacement array, xyz'
			dimslice = slice(0,3)
			distsxyz = [[norm(ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice],axis=1)**2 
				for i in range(nframes-d)] 
				for d in range(nframes)]
			checktime()
		else:
			#---pre-compute a master array of all displacements
			#---Nb the apply_along_axis code is brutally slow so upgrade to numpy 1.8 on dirac
			print 'status: precomputing displacement array, xy'
			dimslice = slice(0,2)
			distsxy = [[norm(ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice],axis=1)**2 
				for i in range(nframes-d)] 
				for d in range(nframes)]
			checktime()
			print 'status: precomputing displacement array, z'
			dimslice = slice(2,3)
			distsz = [[norm(ionstraj[ionsel,i+d,dimslice]-ionstraj[ionsel,i,dimslice],axis=1)**2 
				for i in range(nframes-d)] 
				for d in range(nframes)]
			checktime()
			#---loop over zones
			if zonesub == None:
				zonesub = range(nzones)
			mastermsd_zones = []
			for z in zonesub:
				print 'status: zone = '+str(z)
				inzone = (roundz == z).T
				inzonesliced = [[mean(inzone[i:i+d+1],axis=0) for i in range(nframes-d)] for d in range(nframes)]
				bigmask = [numpy.ma.masked_greater_equal(inzonesliced[i],occupancy) for i in range(nframes)]
				mastermsd = [(array((1*ma.mask_rows(bigmask[i]).mask+1*ma.mask_cols(bigmask[i]).mask)) 
					if shape(bigmask[i].mask) != () else zeros(bigmask[i].shape)) for i in range(nframes)]
				mastermsd_zones.append(mastermsd)
				checktime()
				#---Nb memory flat after a single zone if you delete as follows
				del inzone,inzonesliced,bigmask,mastermsd
			#---Nb currently disabled due to excessive sizes
			if 0:
				#---record times
				times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
				#---package the results
				result_data = MembraneData('ionskate')
				#---Nb data are packaged as type (distsxy,distsz,mastermsdzones)
				#---Nb the mastermsd_zones type goes by zone, ion, delta-t, start frame
				result_data.data = [mastermsd_zones,distsxy,distsz]
				result_data.addnote(['times',times])
				result_data.addnote(['occupancy',occupancy])
				result_data.addnote(['zones',zones])
				result_data.addnote(['ionsel',ionsel])
				#---pickle the results
				pickledump(result_data,resultpkl,directory=pickles)	

#---Nb when selecting subsets, make sure the automatically-generated center gives a symmetric distribution
if 0: plt.grid(True);plt.hist(array(roundz).flatten(),range=(0,30),bins=30);plt.show()

#---POSTPROCESS
#-------------------------------------------------------------------------------------------------------------

#---plotting routine
makeplots = ['panel','calc_diffusion','diffusion_summary','all_raw_msds','mode_diffusion',
	'all_raw_msds_xyz'][1:]

#---postprocess and plot
if 'postproc' in routine:
	#---Nb disabled this to skip writing step because files too big
	if 'skates' not in globals() and 0:
		skates = []
		for aname in analysis_names:
			for i in analysis_descriptors[aname]: vars()[i] = (analysis_descriptors[aname])[i]
			grofile,trajfile = trajectory_lookup(analysis_descriptors,aname,globals())
			resultpkl = 'pkl.ionskate.'+specname_pickle(sysname,trajfile[0])+'.pkl'
			skate = unpickle(pickles+resultpkl)
			skates.append(skate)
		nframes = 500
		times = array(skate.getnote('times'))
		distsxy = skates[0].get(['type','distsxy'])
		distsz = skates[0].get(['type','distsz'])
	#---you have to tune upto according to the time you want and it is modulated by occupancy
	times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)][:upto]
	if 'panel' in makeplots:
		print 'status: making panel plots'
		#---color management
		cslice = int(round(center/bwid))
		if zonesub == None:
			cslice = int(round(center/bwid))
			zone_to_slice = [(i-cslice)+nzones*(-1 if (i-cslice)>(nzones/2) else 0) for i in range(nzones)]
			#---if no subsets, just assume 1 nm bin widths
			plot_slices = [-6,-5,-4,-3,-2,2,3,4,5,6]
		else:
			zone_to_slice = [float(i-cslice)*bwid/10. for i in zonesub]
			plot_slices = [-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11]
			plot_zones = range(len(zonesub))
		redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones for i in range((nzones+1))]]
		bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones for i in range((nzones+1))]]
		clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
		#---figure
		fig = plt.figure()
		gs = gridspec.GridSpec(2,nzones)
		for z in plot_zones:
			print 'zone = '+str(z)
			clr = clrsd[plot_zones.index(z)]
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[0,z])
			allcurves = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsxy[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clr,
					markeredgecolor=clrs[z%len(clrs)],alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_yticklabels([])
			ax.set_xticklabels([])
			ax.set_xlim((10**2,10**3))
			ax.set_ylim((10**-1,10**3))
		for z in plot_zones:
			print 'zone = '+str(z)
			clr = clrsd[plot_zones.index(z)]
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[1,z])
			allcurves = array([mean(ma.masked_values((1*(z0[i]==2)*distsz[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clr,alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_yticklabels([])
			ax.set_xticklabels([])
			ax.set_xlim((10**2,10**3))
			ax.set_ylim((10**-3,10**2))
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-panels-occ'+\
			str(occupancy)+'.png',dpi=300,bbox_inches='tight')
		plt.close(fig)
	#---calculate diffusion rates
	if 'calc_diffusion' in makeplots or ('diffusion_summary' in makeplots and 'diffusexy' not in globals()):
		print 'status: calculating diffusion rates'
		#---fit to get diffusion rates
		diffusexy = []
		diffusez = []
		for z in range(nzones):
			print z
			z0 = mastermsd_zones[z]
			allcurvesxy = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsxy[i]).T,0.).data,axis=1) 
				for i in range(nframes)]).T
			allcurvesz = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsz[i]).T,0.).data,axis=1) 
				for i in range(nframes)]).T
			dconstsxy = []
			dconstsz = []
			if not all(allcurvesxy==0) and not all(allcurvesz==0):
				for c in range(len(allcurvesxy)):
					pair_diffuse = []					
					for curv in [allcurvesxy[c],allcurvesz[c]]:
						curvfilt = array([times,curv[:upto]]).T[curv[:upto]!=0.]
						fitlims = curvfilt[:,0] < timelimit
						if fitlims != []:
							[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
							pair_diffuse.append(bz)
						else: pair_diffuse.append(-1)
					if -1 not in pair_diffuse:
						dconstsxy.append(pair_diffuse[0]/3./2)
						dconstsz.append(pair_diffuse[1]/3./1)
			diffusexy.append(dconstsxy)
			diffusez.append(dconstsz)
	#---summary plot for comparing lateral and normal diffusion
	if 'diffusion_summary' in makeplots:
		print 'status: plotting diffusion rates'
		#---original code before adding zonesub
		if 0:
			#---color management
			cslice = int(round(center/bwid))
			zone_to_slice = [i-cslice for i in range(nzones)]
			zone_to_slice = [(i-cslice)+nzones*(-1 if (i-cslice)>(nzones/2) else 0) for i in range(nzones)]
			plot_slices = [-5,-4,-3,-2,2,3,4,5]
			plot_zones = [zone_to_slice.index(i) for i in plot_slices]
			redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
			bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
			clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
		#---color management
		cslice = int(round(center/bwid))
		if zonesub == None:
			cslice = int(round(center/bwid))
			zone_to_slice = [(i-cslice)+nzones*(-1 if (i-cslice)>(nzones/2) else 0) for i in range(nzones)]
			#---if no subsets, just assume 1 nm bin widths
			plot_slices = [-6,-5,-4,-3,-2,2,3,4,5,6]
			plot_zones = [zone_to_slice.index(i) for i in plot_slices]
			redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones for i in range((nzones+1))]]
			bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones for i in range((nzones+1))]]
			clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
			close_order = range(len(plot_zones))
		else:
			zone_to_slice = [float(i-cslice)*bwid/10. for i in zonesub]
			plot_slices = [-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10]
			plot_zones = range(len(zonesub))
			redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones/1.2 for i in range(nzones-1)]]
			bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones/1.2 for i in range(nzones-1)]]
			clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
			close_order = [nzones/2-1+i*(i%2==0)/2-i*(i%2==1)/2 for i in [int(j) 
				for j in range(1,nzones+1)]][::-1]
		#---generate plot		
		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 40
		allalpha = 1.
		edgeprop = 2.
		#xlims = (10**-4,max([array(i)[array(i)>0].max() for i in diffusexy if i != []])*edgeprop)
		#ylims = (10**-4,max([array(i)[array(i)>0].max() for i in diffusez if i != []])*edgeprop)
		xlims = (10**-4,10**0)
		ylims = (10**-4,10**0)
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in [plot_zones[k] for k in close_order]:
			print z
			clr = clrsd[plot_zones.index(z)]
			dat = array([diffusexy[z],diffusez[z]]).T
			if shape(dat) == (2,): dat = array([dat])
			print max(dat[:,1])
			#---toss any negative slopes
			ax.plot(dat[:,0],dat[:,1],'o',c=clr,markeredgecolor=clr,
				alpha=allalpha,label=str(zone_to_slice[z]))
			#---x axis histogram
			counts,edges = histogram(log10(dat[:,0]),bins=nbins,range=log10(xlims))
			bins = 10**((edges[:-1]+edges[1:])/2.)
			print max(bins)
			axHistx.plot(bins,counts,'-',lw=2.,c=clr,
				markeredgecolor=clr,alpha=allalpha)
			#---y axis histogram
			bins,counts = histogram(log10(dat[:,1]))
			counts,edges = histogram(log10(dat[:,1]),bins=nbins,range=log10(ylims))
			bins = 10**((edges[:-1]+edges[1:])/2.)
			print max(bins)
			axHisty.plot(counts,bins,'-',lw=2.,c=clr,
				markeredgecolor=clr,alpha=allalpha)
		ax.set_xlim(xlims)
		ax.set_ylim(ylims)
		legend_h,legend_l = ax.get_legend_handles_labels()
		ax.legend([legend_h[close_order.index(i)] for i in range(nzones)],[legend_l[close_order.index(i)] 
			for i in range(nzones)],loc='center left',bbox_to_anchor=(1, 0.5),
			bbox_transform=axHisty.transAxes)
		ax.set_yscale('log')
		ax.set_xscale('log')		
		ax.grid(True)
		axHistx.grid(True)
		axHisty.grid(True)
		axHistx.set_xscale('log')
		axHisty.set_yscale('log')		
		ax.set_ylabel(r'$D_Z$',fontsize=fsaxlabel)
		ax.set_xlabel(r'$D_{XY}$',fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		axHistx.set_yticklabels([])
		axHisty.set_xticklabels([])
		plt.setp(axHistx.get_xticklabels()+axHisty.get_yticklabels(),visible=False)
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-diffusions-zdecomp-occ'+\
			str(occupancy)+'.png',dpi=300,bbox_inches='tight')
		plt.show()
	
	#---preferred zone calculation
	if 'all_raw_msds' in makeplots or 'mode_diffusion' in makeplots:
		zonesabs = list(arange(center,0-bwid,-bwid))[::-1][:-1]+list(arange(center,vecs[2]+bwid,bwid))
		zonecenters = array(zonesabs[:-1]+zonesabs[1:])/2.
		rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+\
			(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)]
		rounded = [[int(round(i)) for i in (rewrapz[j]-center)/bwid] for j in range(nions)]
		preferred_zones = stats.mode(rounded,axis=1)[0][:,0]
		rang = int(ptp(preferred_zones))
		pzs = [int(i-preferred_zones.min()) for i in preferred_zones]
		times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
	
	#---plot all the unfiltered MSD curves
	if 'all_raw_msds' in makeplots:
		clrs = [brewer2mpl.get_map('paired','qualitative',12).mpl_colors[i] for i in range(12)]
		fig = plt.figure()
		ax = plt.subplot(121)
		allrawcurves = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
		for c in range(len(allrawcurves)):
			curv = allrawcurves[c]
			ax.plot(times,curv,c=clrs[pzs[c]%len(clrs)])
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim((10**-1,10**6))
		ax = plt.subplot(122)
		allrawcurves = array([mean(distsz[i],axis=0) for i in range(nframes)]).T
		for c in range(len(allrawcurves)):
			curv = allrawcurves[c]
			ax.plot(times,curv,c=clrs[pzs[c]%len(clrs)])
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim((10**-1,10**6))
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-msd-zdecomp-modes.png',
			dpi=300,bbox_inches='tight')
		plt.clf()
		plt.close()
		
	#---plot all the unfiltered MSD curves
	if 'all_raw_msds_xyz' in makeplots:
		clrs = [brewer2mpl.get_map('paired','qualitative',12).mpl_colors[i] for i in range(12)]
		fig = plt.figure()
		ax = plt.subplot(111)
		allrawcurves = array([mean(distsxyz[i],axis=0) for i in range(nframes)]).T
		for c in range(len(allrawcurves)):
			curv = allrawcurves[c]
			ax.plot(times,curv[0:upto],c=clrs[c%len(clrs)]) # It should not be necessary to fix this to upto.
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim((10**-1,10**6))
		plt.show()
#		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-msd-xyz.png',
#			dpi=300,bbox_inches='tight')
		plt.clf()
		plt.close()
		
	#---plot the diffusion summary according to the mode z-slice an ion resides in
	if 'mode_diffusion' in makeplots:
		if 'diffusexy_raw' not in globals() or diffusexy_raw == []:
			print 'status: calculating raw diffusion rates'
			#---fit to get diffusion rates
			allcurvesxy = array([mean(distsxy[i],axis=0) for i in range(nframes)]).T
			allcurvesz = array([mean(distsz[i],axis=0) for i in range(nframes)]).T
			diffusexy_raw = []
			for curv in allcurvesxy[::1]:
				curvfilt = array([times,curv]).T
				fitlims = curvfilt[:,0] < timelimit
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				diffusexy_raw.append(bz/3/2)
			diffusez_raw = []
			for curv in allcurvesz[::1]:
				curvfilt = array([times,curv]).T
				fitlims = curvfilt[:,0] < timelimit
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				diffusez_raw.append(bz/3/1)
			diffusexy_raw = array(diffusexy_raw)
			diffusez_raw = array(diffusez_raw)
		print 'status: plotting raw diffusion rates'
		#---color management
		cslice = int(round(center/bwid))
		zone_to_slice = [i-cslice for i in range(nzones)]
		zone_to_slice = [(i-cslice)+nzones*(-1 if (i-cslice)>(nzones/2) else 0) for i in range(nzones)]
		plot_slices = [-5,-4,-3,-2,2,3,4,5]
		plot_zones = [zone_to_slice.index(i) for i in plot_slices]
		redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
		bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
		clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
		#---p
		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 40
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in plot_zones:
			clr = clrsd[plot_zones.index(z)]
			dat = array([diffusexy_raw[array(pzs)==z],diffusez_raw[array(pzs)==z]]).T	
			if shape(dat) == (2,): dat = array([dat])
			#---toss any negative slopes
			dat = array([i for i in dat if i[0] > 0 and i[1] > 0])
			ax.plot(dat[:,0],dat[:,1],'o',c=clr,markeredgecolor=clr,
				alpha=allalpha,label=str(zone_to_slice[z]))
			#---x axis histogram
			counts,edges = histogram(log10(dat[:,0]),bins=nbins,range=log10(xlims))
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHistx.plot(bins,counts,'-',lw=2.,c=clr,
				markeredgecolor=clr,alpha=allalpha)
			#---y axis histogram
			bins,counts = histogram(log10(dat[:,1]))
			counts,edges = histogram(log10(dat[:,1]),bins=nbins,range=log10(ylims))
			bins = 10**((edges[:-1]+edges[1:])/2.)
			axHisty.plot(counts,bins,'-',lw=2.,c=clr,
				markeredgecolor=clr,alpha=allalpha)	
		edgeprop = 2.
		xlims = (min(diffusexy_raw)/edgeprop,max(diffusexy_raw)*edgeprop)
		ylims = (min(diffusez_raw)/edgeprop,max(diffusez_raw)*edgeprop)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),bbox_transform=axHisty.transAxes)
		ax.set_ylim(ylims)
		ax.set_xlim(xlims)
		axHisty.set_ylim(ylims)
		axHistx.set_xlim(xlims)
		ax.set_yscale('log')
		ax.set_xscale('log')		
		ax.grid(True)
		axHistx.grid(True)
		axHisty.grid(True)
		axHistx.set_xscale('log')
		axHisty.set_yscale('log')		
		ax.set_ylabel(r'$D_Z\:(units)$',fontsize=fsaxlabel)
		ax.set_xlabel(r'$D_{XY}\:(units)$',fontsize=fsaxlabel)
		plt.setp(ax.get_yticklabels(),fontsize=fsaxlabel)
		plt.setp(ax.get_xticklabels(),fontsize=fsaxlabel)
		axHistx.set_yticklabels([])
		axHisty.set_xticklabels([])X
		plt.setp(axHistx.get_xticklabels()+axHisty.get_yticklabels(),visible=False)
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-diffusions-zdecomp-modes.png',
			dpi=300,bbox_inches='tight')
		plt.show()
		
	#---the following code 
	if 'check_timescales' in makeplots:
		times = [mean([clock[i+d]-clock[i] for i in range(nframes-d)]) for d in range(nframes)]
		#---select zone/slice -5
		z = 0
		#---grab the slow indices for plotting
		slow_inds = [i for i in range(len(diffusexy[0])) if diffusexy[0][i] < 10**-2]
		fast_inds = [i for i in range(len(diffusexy[0])) if diffusexy[0][i] > 10**-2]
		z0 = mastermsd_zones[z]
		#---Nb the occupancy diffusion calculations were missing the equals two criteria
		#allcurves = array([mean(ma.masked_values((z0[i]*distsxy[i]).T,0.).data,axis=1) 
		#	for i in range(nframes)]).T
		allcurves = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsxy[i]).T,0.).data,axis=1) 
			for i in range(500)]).T
		fig = plt.figure()
		ax = plt.subplot(111)
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		for i in slow_inds:
			ax.plot(times,allcurves[i],'b')
		for i in fast_inds:
			ax.plot(times,allcurves[i],'r')
		plt.show()
		


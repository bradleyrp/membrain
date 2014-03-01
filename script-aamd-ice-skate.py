#!/usr/bin/python

if 'mset' not in globals():
	interact = True
	from membrainrunner import *
	execfile('locations.py')

from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---Settings
#-------------------------------------------------------------------------------------------------------------

#---possible analyses
analysis_descriptors = {
	'v511-30000-80000-100':
		{'sysname':'membrane-v511',
		'sysname_lookup':'membrane-v511-ions',
		'trajsel':'s6-kraken-md.part0009.30000-80000-100.ions.xtc',
		'structure_pkl':
			'pkl.structures.membrane-v511.a2-surfacer.s6-kraken-md.part0009.30000-80000-100.pkl',
		'ionname':'Cal'}}
analysis_names = ['v511-30000-80000-100']
routine = ['compute','postproc'][1:]

#---method
upto = 10
occupancy = 0.9
bwid = 10

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---decomposing the diffusion by calculating it only from movements that are mostly in a particular z-slice
if 'compute' in routine or 'mastermsd_zones' not in globals():
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
		ionspos = array(ionspos)[:-1]
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
		thick = mean(mset_surf.surf_thick)
		zonesabs = list(arange(center,0-bwid,-bwid))[::-1][:-1]+list(arange(center,vecs[2]+bwid,bwid))
		zonecenters = array(zonesabs[:-1]+zonesabs[1:])/2.
		rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+\
			(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)]
		roundz = [[int(round(i)) for i in (rewrapz[j]-center)/bwid] for j in range(nions)]
		roundz = roundz - array(roundz).min()
		nzones = ptp(roundz)
		#---select subset of ions if desired
		ionsel = slice(0,nions)
		nions = len(range(nions)[ionsel])
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
		mastermsd_zones = []
		for z in range(nzones):
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

#---POSTPROCESS
#-------------------------------------------------------------------------------------------------------------

#---plotting routine
makeplots = ['panel','diffusion_summary','all_raw_msds','mode_diffusion'][1:2]

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
		fig = plt.figure()
		gs = gridspec.GridSpec(2,nzones)
		for z in range(nzones):
			print 'zone = '+str(z)
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[0,z])
			allcurves = array([mean(ma.masked_values(array(1*(z0[i]==2)*distsxy[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clrs[z%len(clrs)],
					markeredgecolor=clrs[z%len(clrs)],alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
		for z in range(nzones):
			print 'zone = '+str(z)
			z0 = mastermsd_zones[z]
			ax = fig.add_subplot(gs[1,z])
			allcurves = array([mean(ma.masked_values((1*(z0[i]==2)*distsz[i]).T,0.).data,axis=1) 
				for i in range(upto)]).T
			for curv in allcurves[::1]:
				valids = curv != 0.
				curvfilt = array([times,curv]).T[valids]
				ax.plot(curvfilt[:,0],curvfilt[:,1],'-',c=clrs[z%len(clrs)],alpha=0.7)
			ax.set_xscale('log')
			ax.set_yscale('log')
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-panels-occ'+\
			str(occupancy)+'.png',dpi=300,bbox_inches='tight')
		plt.clf()
		plt.close()
	#---calculate diffusion rates
	if 'calc_diffusion' in makeplots or ('diffusion_summary' in makeplots and 'diffusexy' not in globals()):
		print 'status: calculating diffusion rates'
		#---fit to get diffusion rates
		diffusexy = []
		for z in range(nzones):
			print z
			z0 = mastermsd_zones[z]
			allcurves = array([mean(ma.masked_values((z0[i]*distsxy[i]).T,0.).data,axis=1) 
				for i in range(nframes)]).T
			dconsts = []
			if not all(allcurves==0):
				for curv in allcurves[::1]:
					curvfilt = array([times,curv[:upto]]).T[curv[:upto]!=0.]
					fitlims = curvfilt[:,0]<1*10**10
					[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
					dconsts.append(bz/3/2)
			diffusexy.append(dconsts)
		diffusez = []
		for z in range(nzones):
			print z
			z0 = mastermsd_zones[z]
			allcurves = array([mean(ma.masked_values((z0[i]*distsz[i]).T,0.).data,axis=1) 
				for i in range(nframes)]).T
			dconsts = []
			if not all(allcurves==0):
				for curv in allcurves[::1]:
					curvfilt = array([times,curv[:upto]]).T[curv[:upto]!=0.]
					fitlims = curvfilt[:,0]<1*10**3
					[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
					dconsts.append(bz/3/1)
			diffusez.append(dconsts)
	#---summary plot for comparing lateral and normal diffusion
	if 'diffusion_summary' in makeplots:
		print 'status: plotting diffusion rates'
		#---color management
		cslice = int(round(center/bwid))
		zone_to_slice = [i-cslice for i in range(nzones)]
		zone_to_slice = [(i-cslice)+nzones*(-1 if (i-cslice)>(nzones/2) else 0) for i in range(nzones)]
		plot_slices = [-5,-4,-3,-2,2,3,4,5]
		plot_zones = [zone_to_slice.index(i) for i in plot_slices]
		redclrs = ['k']+[mpl.cm.RdBu_r(j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
		bluclrs = ['k']+[mpl.cm.RdBu_r(1.-j) for j in [float(i)/nzones for i in range((nzones+1)/2)]]
		clrsd = [(redclrs[-i] if i < 0 else bluclrs[i]) for i in plot_slices]
		#---generate plot		
		fig = plt.figure()
		ax = plt.subplot(111)
		nbins = 20
		allalpha = 1.
		edgeprop = 2.
		xlims = (min(diffusexy[diffusexy>0])/edgeprop,max(diffusexy[diffusexy>0])*edgeprop)
		ylims = (min(diffusez[diffusexy>0])/edgeprop,max(diffusez[diffusexy>0])*edgeprop)
		divider = make_axes_locatable(ax)
		axHistx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
		axHisty = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
		for z in plot_zones:
			clr = clrsd[plot_zones.index(z)]
			dat = array([diffusexy[z],diffusez[z]]).T
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
		nonzeroxy = [i for i in flatten(array(diffusexy)) if i > 0.]
		nonzeroz = [i for i in flatten(array(diffusez)) if i > 0.]
		xlims = (min(diffusexy[diffusexy>0])/edgeprop,max(diffusexy[diffusexy>0])*edgeprop)
		ylims = (min(diffusez[diffusexy>0])/edgeprop,max(diffusez[diffusexy>0])*edgeprop)
		ax.set_ylim(ylims)
		ax.set_xlim(xlims)
		axHisty.set_ylim(ylims)
		axHistx.set_xlim(xlims)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),bbox_transform=axHisty.transAxes)
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
				fitlims = curvfilt[:,0]<1*10**3
				[bz,az] = numpy.polyfit(curvfilt[fitlims,0],curvfilt[fitlims,1],1)
				diffusexy_raw.append(bz/3/2)
			diffusez_raw = []
			for curv in allcurvesz[::1]:
				curvfilt = array([times,curv]).T
				fitlims = curvfilt[:,0]<1*10**3
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
		axHisty.set_xticklabels([])
		plt.setp(axHistx.get_xticklabels()+axHisty.get_yticklabels(),visible=False)
		plt.savefig(pickles+'fig-'+sysname.split('-')[-1]+'-iceskate-diffusions-zdecomp-modes.png',
			dpi=300,bbox_inches='tight')
		plt.show()
		


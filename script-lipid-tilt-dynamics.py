#!/usr/bin/python -i

if 0:

	from membrainrunner import *

	#---Settings
	#-------------------------------------------------------------------------------------------------------------

	#---Analysis parameters
	skip = None
	framecount = None
	location = ''
	execfile('locations.py')

	#---Selections
	director_cgmd = ['name PO4','name C4A','name C4B']
	selector_cgmd = 'name PO4'
	cgmd_protein = 'name BB'

	#---Analysis plan
	analysis_plan = slice(0,None)
	analysis_descriptors = [
		(['membrane-v700'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
		(['membrane-v701'],director_cgmd,selector_cgmd,cgmd_protein,slice(-1,None)),
		(['membrane-v550'],director_cgmd,selector_cgmd,None,slice(-1,None))]
	
	#---MAIN
	#-------------------------------------------------------------------------------------------------------------

	'''
	outline
	open tilt calculations
	make an mset for the locations
	for a single frame, plot positions and color
	'''

	msetset = []

	starttime = time.time()
	print 'Starting analysis job.'
	for ad in analysis_descriptors[analysis_plan]:
		(tests,director,selector,protein_selection,trajno) = ad
		sysname = tests[0]
		residues = selector
		for t in range(len(tests)):
			testno = t
			print 'Loading a membrane set '+tests[t]+'.'
			for traj in trajectories[systems.index(tests[t])][trajno]:
				mset = MembraneSet()
				gro = structures[systems.index(tests[testno])]
				basename = traj.split('/')[-1][:-4]
				sel_surfacer = sel_aamd_surfacer
				print 'Accessing '+basename+'.'
				mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
					resolution='cgmd')
				#---Average structure calculation
				mset.identify_monolayers(director,startframeno=0)
				end = None
				start = None
				if framecount == None:
					if end == None: end = mset.nframes
					if start == None: start = 0
					if skip == None: skip = 1
				else:
					start = 0
					end = mset.nframes
					skip = int(float(mset.nframes)/framecount)
					skip = 1 if skip < 1 else skip
				mset_tilts = unpickle(pickles+'pkl.tilt.'+sysname[9:14]+'.'+basename+'.pkl')
				mset.store.append(mset_tilts.getdata('tilts'))
				del mset_tilts
				msetset.append(mset)
	#print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'
if 1:
	titles = ['parallel','antiparallel','control']
	plt.rc('font', family='sans-serif')
	plot_video_filebase = 'vid-lipid-tilts-v700.v701.v550'
	angle_thresh = 130.
	framenums = [i[0] for i in msetset[0].getdata('tilts').label]
	#---note about frames here
	for fr in framenums:
		print 'processing frame '+str(fr)
		fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(15,5))
		for sys in range(3):
			msetset[sys].gotoframe(fr)
			topxyz = array([mean(msetset[sys].universe.residues[i].selectAtoms(selector).coordinates(),axis=0) 
				for i in msetset[sys].monolayer_residues[0]])
			vecs = msetset[sys].vec(fr)
			axes[sys].set(aspect=1)
			axes[sys].set_xlim((0,vecs[0]))
			axes[sys].set_ylim((0,vecs[0]))
			taila = [180./pi*k if (not isnan(k)) else 180. 
				for k in array(flatten(msetset[sys].getdata('tilts').get(['monolayer',0,'tail',0,'frame',framenums.index(fr)])))]
			tailb = [180./pi*k if (not isnan(k)) else 180. 
				for k in array(flatten(msetset[sys].getdata('tilts').get(['monolayer',0,'tail',1,'frame',framenums.index(fr)])))]
			tilts = mean([array(taila),array(tailb)],axis=0)
			patches = []
			patches_special = []
			for pt in topxyz:
				patches.append(plt.Circle([pt[0],pt[1]],alpha=0.5,radius=5.0))
			patches_special = [plt.Circle([topxyz[p][0],topxyz[p][1]],alpha=0.5,radius=8.0) 
				for p in range(len(topxyz)) if tilts[p] < angle_thresh]
			patchcoll = mpl.collections.PatchCollection(patches,cmap=mpl.cm.jet,linewidths=0.)
			patchcoll.set_array(np.array(tilts))
			patchcoll.set_alpha(0.5)
			axes[sys].add_collection(patchcoll)
			patchcoll2 = mpl.collections.PatchCollection(patches_special,cmap=mpl.cm.gray,linewidths=0.)
			patchcoll2.set_alpha(1.0)
			patchcoll2.set_array(np.array([1. for i in range(len(patches_special))]))
			axes[sys].add_collection(patchcoll2)
			if analysis_descriptors[sys][3] != None:
				protpts = msetset[sys].universe.selectAtoms(analysis_descriptors[sys][3]).coordinates()[:,0:2]
				hull = scipy.spatial.ConvexHull(protpts)
				axes[sys].plot(protpts[hull.vertices,0],protpts[hull.vertices,1],'r-',lw=1)
			axes[sys].set_title(titles[sys])
		plt.savefig(pickles+plot_video_filebase+'.fr.'+str('%04d'%framenums.index(fr))+'.png',dpi=500,bbox_inches='tight')
		plt.cla()
		plt.close()
	subprocess.call(['ffmpeg','-i',pickles+'/'+plot_video_filebase+'.fr.%04d.png','-vcodec','mpeg2video','-qscale','0','-filter:v','setpts=2.0*PTS',pickles+'/'+plot_video_filebase+'.mpeg'])
	#os.popen('rm -r -f '+pickles+'/figs-'+sysname+'-tilefilter.snapshots')


#!/usr/bin/python -i

from membrainrunner import *

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---method
skip = None
framecount = 200
framecount = None # Pickle will have same slicing as gmx-time-slice
framecount_force_exact = True
location = ''
execfile('locations.py')

#---selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'
sel_aamd_surfacer = ['name P','(name C2 and not resname CHL1)']

#---analysis plan
analysis_descriptors = [
	(['membrane-v531-gr2d'],'all',director_asymmetric,[-1],
		[['resname POPC and name P','resname POPC and name P','POPC-POPC']],'v531.POPC-POPC'),
#	(['membrane-v532'],'all',director_asymmetric,[-2,-1],
#		[['resname POPC and name P','resname POPC and name P','POPC-POPC']],'v532.POPC-POPC')
	]
analyses = [analysis_descriptors[i] for i in [-1]]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over analysis descriptors
for ad in analyses:
	(tests,selector,director,trajnos,pairs,exptname) = ad
	#---loop over tests within the descriptor
	for testno in range(len(tests)):
		#---loop over specified trajectories
		for traj in [trajectories[systems.index(tests[testno])][i] for i in trajnos]:
			mset = MembraneSet()
			#---load
			gro = structures[systems.index(tests[testno])]
			basename = traj.split('/')[-1][:-4]
			sel_surfacer = sel_aamd_surfacer
			print 'Accessing '+basename+'.'
			mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
				resolution='aamd')
			mset.identify_monolayers(director)
			#---frame selection header
			end = None
			start = None
			if framecount == None:
				if end == None: end = mset.nframes
				if start == None: start = 0
				if skip == None: skip = 1
				framerange = range(start,end,skip)
			else:
				start = 0
				end = mset.nframes
				skip = int(round(float(mset.nframes)/framecount))
				skip = 1 if skip < 1 else skip
				if framecount_force_exact:
					framerange = list([int(round(i)) for i in linspace(0,mset.nframes-1,framecount)])
				else:
					framerange = range(start,end,skip)
			for pair in pairs:
				result_data = MembraneData('gr2d',label=pair[2])
				allcurves = []
				for mononum in range(2):
					for lnum in range(2):
						allselect_lipids = mset.universe.selectAtoms(pair[lnum])
						validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
							set([i-1 for i in allselect_lipids.resids()])))
						mset.selections.append(sum([allselect_lipids.residues[
							list(allselect_lipids.resids()).index(i+1)].selectAtoms(pair[lnum]) 
							for i in validresids]))
					allcurves_by_monolayer = []
					if mset.selections[0] == 0.0 or mset.selections[1] == 0.0:
						print 'Note: missing lipids in this monolayer.'
						allcurves_by_monolayer = [[] for i in range(start,end,skip)]
					else:
						for frameno in framerange:
							mset.gotoframe(frameno)
							vecs = mset.vec(frameno)
							pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
							pts2 = array(mset.get_points(frameno,selection_index=1))[:,0:2]
							points = pts2
							dims=[0,1]
							ans = []
							for p in points:
								for tr in [[i,j] for i in arange(-2,2+1,1) for j in arange(-2,2+1,1) 
									if not (i == 0 and j == 0)]:
									ans.append([p[i]+tr[i]*vecs[i] if i in dims else p[i] for i in range(2)])
							pts2pbc = concatenate((points,array(ans)))
							cutoff = 2*(vecs[0] if vecs[0]<vecs[1] else vecs[1])
							sysarea = pi*cutoff**2
							dmat2 = scipy.spatial.distance.cdist(pts1,pts2pbc)
							binsizeabs = 4
							scanrange = arange(0,int(cutoff),binsizeabs)
							histcomb = []
							for r in range(len(dmat2)):
								row = dmat2[r]
								hist,binedge = numpy.histogram(array(row)[array(row)!=0.],
									range=(0.,int(max(scanrange))),bins=len(scanrange)-1)
								mid = (binedge[1:]+binedge[:-1])/2
								histcomb.append(hist)
							grcurve = sum(histcomb,axis=0)/float(len(histcomb))
							allcurves_by_monolayer.append(grcurve)
						cutoff = 2*min([int(i) for i in np.mean(mset.vecs,axis=0)[0:2]])
					allcurves.append(allcurves_by_monolayer)
					mset.selections = []
				points_counts = [shape(pts1),shape(pts2pbc)]
				pair_selects = pair[0:2]
				pair_name = pair[2]
				#---load up the membranedata object with calculation details
				savelist = ['cutoff','sysarea','points_counts','binsizeabs','pair_selects','pair_name']
				for item in savelist:
					result_data.addnote([item,globals()[item]])
				result_data.data = [[allcurves[0][i],allcurves[1][i]] for i in range(len(framerange))]
				result_data.label = framerange
				mset.store.append(result_data)
				del result_data
			pickledump(mset,'pkl.gr2d.'+exptname+'.'+basename.strip('md.')+'.pkl',directory=pickles)
			del mset
			
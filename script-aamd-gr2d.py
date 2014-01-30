#!/usr/bin/python -i

from membrainrunner import *

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

#---settings
skip = None
framecount = None
location = ''
execfile('locations.py')

#---selections
director_symmetric = ['name P','name C218','name C318']
director_asymmetric = ['(name P and not resname CHL1) or (name C3 and resname CHL1)',
		'(name C218 and not resname CHL1) or (name C25 and resname CHL1)']
selector = '(name P and not resname CHL1) or (name C3 and resname CHL1)'

#---analysis plan
analysis_descriptors = [
	(['membrane-v510'],'all',director_symmetric,-1,
		[['resname DOPC and name P','name CL','DOPC P-to-CL'],
		['resname DOPC and name P','name MG','DOPC P-to-MG']]),
	(['membrane-v530'],'all',director_asymmetric,-1,
		[['resname DOPS and name P','resname PI2P and name P','DOPS-PIP2']]),
	(['membrane-v509'],'all',director_symmetric,-1,
		[['resname PI2P and name P','resname PI2P and name P','DOPS-DOPS']]),
	(['membrane-v530'],'all',director_asymmetric,-1,
		[['resname PI2P and name P','resname PI2P and name P','PIP2-PIP2']])]
analyses = [analysis_descriptors[i] for i in [-1]]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---loop over analysis descriptors
for ad in analyses:
	(tests,selector,director,trajno,residues,pairs,extraname) = ad
	#---loop over tests within the descriptor
	for testno in range(len(tests)):
		#---loop over specified trajectories
		for traj in trajectories[systems.index(tests[testno])][trajno]:
			(tests,selector,director,trajno,pairs) = analysis_descriptors[-1]
			mset = MembraneSet()
			#---frame selection header
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
			#---load
			gro = structures[systems.index(tests[testno])]
			basename = traj.split('/')[-1][:-4]
			sel_surfacer = sel_aamd_surfacer
			print 'Accessing '+basename+'.'
			mset.load_trajectory((basedir+'/'+gro,basedir+'/'+traj),
				resolution='aamd')
			mset.identify_monolayers(director)
			for pair in pairs:
				for lnum in range(2):
					allselect_lipids = mset.universe.selectAtoms(pair[lnum])
					for mononum in range(2):
						validresids = list(set.intersection(set(mset.monolayer_residues[mononum]),
							set([i-1 for i in allselect_lipids.resids()])))
						mset.selections.append(sum([allselect_lipids.residues[
							list(allselect_lipids.resids()).index(i+1)].selectAtoms(pair[lnum]) 
							for i in validresids]))
				allcurvs = []
				for frameno in range(start,end,skip):
					vecs = mset.vec(frameno)
					pts1 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
					pts2 = array(mset.get_points(frameno,selection_index=0))[:,0:2]
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
					hist,binedge = numpy.histogram(array(dmat2[0])[array(dmat2[0])!=0.],
						range=(0.,int(cutoff)),bins=int(cutoff/binsizeabs))
					mid = (binedge[1:]+binedge[:-1])/2
					binwidth = mid[1]-mid[0]
					areas = [pi*binwidth*mid[i]*2 for i in range(len(binedge)-1)]
					areas = [pi*((binedge[i+1])**2-(binedge[i])**2) for i in range(len(binedge)-1)]
					histcomb = []
					for r in range(len(dmat2)):
						row = dmat2[r]
						hist,binedge = numpy.histogram(array(row)[array(row)!=0.],
							range=(0.,int(cutoff)),bins=int(cutoff/binsizeabs))
						histcomb.append(hist)
					grcurve = sum(histcomb,axis=0)/float(len(histcomb))
					allcurvs.append(grcurve)
				mset.selections = []
	
			#---write output ??????????????
			result_data = MembraneData('gr2d',label=)
			mset.store.append(result_data)
			del result_data
			
'''
#---plots
mid = (binedge[1:]+binedge[:-1])/2
binwidth = mid[1]-mid[0]
plt.plot(mid,mean(allcurvs,axis=0)/areas/(len(dmat2)/(vecs[0]*vecs[1])),'o-',c='r',label='reg')
plt.xlim((0,100))
plt.ylim((0,2))
plt.show()
'''


#!/usr/bin/python -i

from membrainrunner import *

location = 'light'
execfile('locations.py')

execfile('plotter.py')
import scipy.interpolate
import scipy.integrate
import os

#---PARAMETERS
#-------------------------------------------------------------------------------------------------------------

#---Parameters, initial (best)
nnnum = 5 			#---number of nearest neighbors
numgridpts = 15 	#---how many grid points, per direction, to sample (20 is optimal)
distance_factor = 2	#---exponent on distance scaling in pseudo-RBF
span = 3			#---how many voxels to move

#---Parameters, fixed
voxelsize=1.0		#---voxel size from the rerun
flatref = False		#---use a flat reference midplane instead of the rounded average structure
view3d = False		#---view the map in 3D using mayavi
exag = 0 			#---exaggerates the peaks on the plot if you use the view3d
check_surf_mean = 0	#---whether to print the mean surface for 
protlenshow = 22 	#---Include only the first few protein points in the centroid calculation (o/w None)
plotunsmoothed = 0	#---Plot the raw data
brokeversion = 0	#---A weird effect with small spans and the wrong input to rezip grid
plotvisible = 0		#---Show the plots if you are doing the comparison
avgsurfround = 0	#---Whether you round the average surface to the voxel grid

#---Parameters, sweep
scan_span = [6,5,4,3,2,1]
scan_nnum = [4,16,32,64,1]
scan_numgridpts = [64,32,20,16]
scan_distance_factor = [2,1]

#---Parameters, sweep
scan_span = [4,3]
scan_nnum = [4,16,32,64,1]
scan_numgridpts = [64,32,20,16]
scan_distance_factor = [1]

#---Parameters, sweep
scan_span = [3,]
scan_nnum = [16]
scan_numgridpts = [64]
scan_distance_factor = [1]

#---Analysis plan
analysis_plan = slice(None)
analysis_descriptors = [
	['v614.part0002','pkl.structures.membrane-v614-stress.md.part0002.rerun.pkl',
		'localpressure.v614.part0002.3Dpp.dat',4],
	['v612.part0003	','pkl.structures.membrane-v612-stress.md.part0003.rerun.pkl',
		'localpressure.v612.part0003.3Dpp.dat',2],
	['v550.part0008','pkl.structures.membrane-v550.md.part0008.rerun.pkl',
		'localpressure.v550.part0008.3Dpp.dat',0],
	['v550.part0008','pkl.structures.membrane-v550.md.part0008.rerun.pkl',
		'localpressure.v550.part0008.3Dpp.dat',0]]
tests = [[3,32,16,1],
	[3,32,64,1],
	[3,32,20,1],
	[3,4,20,1]]
tests = [tests[-2]]
#---Selected parameters which look good
tests = [
	[3,2,16,1],
	[4,32,64,1],
	[3,4,64,1],
	[3,32,64,1],
	[3,64,16,1]]
datdir3dpp = 'localpressure.v614.framewise'
datdir3dpp = 'localpressure.v550.framewise'

#---Parameters, sweep on the best one from above
scan_span = [4]
#scan_nnum = [2,4,5,8,10,12,16,18,20,24,28,32,40,50,64,128,196,250,400]
scan_nnum = [2,4,16,32,64,128]
scan_numgridpts = [64]
scan_distance_factor = [1]
ordering = [scan_span,scan_nnum,scan_numgridpts,scan_distance_factor]
tests = [[i,j,k,l] for i in ordering[0] for j in ordering[1] for k in ordering[2] for l in ordering[3]]

#---Storing the results
storedir = pickles+'localpressure.results.2013.11.27.1800.test'
logfile = pickles+'localpressure.results.2013.11.27.1800.test'+'/'+\
	'log.localpressure.results.2013.11.27.1800'

#---Type of looping or parameter sweeps to do
batch_parameter_sweep = False
batch_comparison = False
batch_parameter_sweep_framewise = True

mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams.update({'font.style':'sans-serif'})
mpl.rcParams.update({'font.size': 16})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	
#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def calculate_stressmap(span,nnnum,numgridpts,distance_factor,imagefile=None,runplots=True,logfile=None,
	plotunsmoothed=False,brokeversion=False):
	'''Curvature calculation via stress tensor, post-post processing.'''
	for pt in dat3dpp:
		if ((abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span) and 
			(abs(pt[2]-avginterpsurf[pt[0]][pt[1]]) < span)):
			rawresults[int(pt[0])][int(pt[1])].append((1./2*(pt[3]+pt[7])-pt[11])*
				(pt[2]-avginterpsurf[pt[0]][pt[1]])*blocksize[2]/10*voxelsize)
	results = array([[scipy.integrate.simps(rawresults[i][j]) for j in range(griddims[1]+1)] 
		for i in range(griddims[0]+1)])
	print shape(results)
	#---Pseudo-Radial basis function computation with K-nearest neighbors
	pts = mset.wrappbc(mset.unzipgrid(results/10**(-exag),vecs=mset.vec(0)),
		vecs=mset.vec(0),mode='grow',growsize=1.0)
	print shape(pts)
	#meshplot(pts,vecs=vecs,show='surf')
	#raw_input('...')
	newpts = array([[i,j] for i in linspace(-vecs[0],2*vecs[0],3*numgridpts) 
		for j in linspace(-vecs[1],2*vecs[1],3*numgridpts)])
	print shape(newpts)
	#---Distance matrix
	tree = scipy.spatial.cKDTree(pts[:,0:2])
	#---Unpack the points and perform the pseudo-RBF
	smoothed = []
	for pt in newpts:
		tmp = tree.query(pt,nnnum)
		#print tmp
		#print shape(tmp[0])
		smoothed.append([pt[0],pt[1],mean([pts[tmp[1][i],2]*1/((tmp[0][i])*sum(tmp[0]))**distance_factor 
			for i in range(len(tmp[0]))])])
	print 'shape(smoothed) = '+str(shape(smoothed))
	#nearset = [i for i in smoothed[1:] if np.linalg.norm(i[0:2]-mean(mset.protein[0],axis=0)[0:2])<300.]
	#---Drop infinite values and the edges, which are noisy
	nearset = [i for i in smoothed if np.isinf(i[2]) == False 
		and (i[0] != 0.0 and i[1] != 0.0 and i[0] != vecs[0] and i[1] != vecs[1])
		and (i[0] > 0. and i[0] < vecs[0] and i[1] > 0. and i[1] < vecs[1])]
	print 'shape(nearset) = '+str(shape(nearset))
	if view3d:
		meshplot(nearset,show='surf')
	if view3d and nprots > 0:
		meshpoints(mset.protein[0]+[-vecs[0]*0,0,-mean(mset.surf_position)-25+50],scale_factor=10,
			color=(0,0,0),opacity=0.2)
	#---Report maximum and mean curvatures assuming kappa = 20 kBT
	peakmax = (max(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	peakmin = (min(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	peakmean = (mean(array(nearset[1:])[:,2])*10**exag)*(10**5*10**(-9*2))/(20*1.38*10**-23*310*10**9)
	print 'peakmax = '+str(peakmax)
	print 'peakmin = '+str(peakmin)
	print 'peakmean = '+str(peakmean)
	logmaxminmean.append([span,nnum,numgridpts,distance_factor,peakmax,peakmin,peakmean])
	global protlenshow
	if nprots > 0:
		protlen = int(shape(mset.protein[0])[0]/nprots)
		if protlenshow == None:
			protlenshow = protlen			
		protein_centers = mean([mean(mset.protein,axis=0)[i*protlen:i*protlen+protlenshow] 
			for i in range(nprots)],axis=1)
	else:
		protein_centers = None
	if plotunsmoothed == True:
		dat = results
	elif brokeversion == True:
		dat = mset.rezipgrid(nearset,vecs=vecs)
	else:
		dat = mset.rezipgrid(nearset,vecs=vecs,grid=[numgridpts,numgridpts])
	return (dat,protein_centers)
	
def plot_stressmap(dat,protein_centers,nprots,numgridpts,imagefile=None,plotvisible=False):
	'''Plot a single stress map.'''
	fig = plt.figure()	
	ax = plt.subplot2grid((1,4),(0,0))
	plt.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft')
	if nprots > 0:
		for pt in protein_centers:
			circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
				int(round(pt[1]/vecs[0]*numgridpts))),radius=0.25,color='k',
				alpha=1.0)
			ax.add_patch(circ)
	if imagefile != None:
		plt.savefig(imagefile)
	if plotvisible:
		plt.show()
	else:
		plt.clf()

#---MAIN
#-------------------------------------------------------------------------------------------------------------
		
#---Compute comparison plots
if batch_comparison:
	starttime = time.time()
	print 'Starting analysis job.'
	#---Sweep parameters and render plots
	for t in range(len(tests)):
		print 'Running postprocessing: stress tensor map scan '+str(tests[t])+'.'
		span,nnum,numgridpts,distance_factor = tests[t]
		result_stack = []
		for ad in analysis_descriptors[analysis_plan]:
			#---Load global variables with calculation specifications used in analysis functions above.
			(systemname,msetfile,picklefile,nprots) = ad
			print 'Running postprocessing: stress tensor map system '+systemname+'.'
			mset = unpickle(pickles+msetfile)
			file3dpp = pickles+picklefile
			#---Make a directory for figures
			if not os.path.isdir(storedir):
				os.mkdir(storedir)
			logmaxminmean = []
			#---Perform global processing steps for this system
			dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
			mset.calculate_average_surface()
			vecs = mean(mset.vecs,axis=0)
			griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
			rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
			xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
				griddims[1]+2)])
			unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
			interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],
				fill_value=0.0)
			blocksize = (vecs/griddims)
			themeanpos = mean(mset.surf_position)
			if flatref == False:
				if avgsurfround:
					avginterpsurf = array([[round((interp(i,j)+themeanpos)/blocksize[2])
						for i in linspace(0,vecs[0],griddims[0]+2)] 
						for j in linspace(0,vecs[1],griddims[1]+2)])
				else:
					avginterpsurf = array([[((interp(i,j)+themeanpos)/blocksize[2])
						for i in linspace(0,vecs[0],griddims[0]+2)] 
						for j in linspace(0,vecs[1],griddims[1]+2)])
			else:
				avginterpsurf = array([[themeanpos/blocksize[2]
					for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
			if ((themeanpos/blocksize[2] - span < 0.0) or 
				(themeanpos/blocksize[2] + span >  max(dat3dpp[:,2]))):
				print 'Warning: your span exceeds your box dimensions'
			#---Done global steps, now proceed to sweep analysis parameters
			print 'Running postprocessing: stress tensor map scan '+str(tests[t])+'.'
			span,nnum,numgridpts,distance_factor = tests[t]
			result_stack.append(calculate_stressmap(span,nnum,numgridpts,distance_factor,
				plotunsmoothed=plotunsmoothed,brokeversion=brokeversion,
				imagefile=storedir+'/stressmap.'+systemname+
				'-span'+str(span)+
				'-nnum'+str(nnum)+
				'-ngp'+str(numgridpts)+
				'-df'+str(distance_factor)+
				'.png'))
		print '... '+str(numgridpts)+' ...'
		fig = plt.figure()	
		gs = mpl.gridspec.GridSpec(4,1,width_ratios=[1,1,2,1],height_ratios=[1])
		plt.rc('font', family='sans-serif')
		extremum = max([max([max(i) for i in result_stack[j][0]]) for j in range(3)])
		ax0 = plt.subplot2grid((1,4),(0,0))
		ax0.set_title(r'$\textbf{{ENTH}\ensuremath{\times}4}$')
		ax0.set_xticklabels([])
		ax0.set_yticklabels([])
		ax0.set_adjustable('box-forced')
		dat = result_stack[0][0]
		protein_centers = result_stack[0][1]
		ax0.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		nprots = 4
		if nprots > 0:
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax0.add_patch(circ)
		ax1 = plt.subplot2grid((1,4),(0,1))
		ax1.set_title(r'$\textbf{{ENTH}\ensuremath{\times}1}$')
		ax1.set_xticklabels([])
		ax1.set_yticklabels([])
		ax1.set_adjustable('box-forced')
		dat = result_stack[1][0]
		protein_centers = result_stack[1][1]
		ax1.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		nprots = 2
		if nprots > 0:
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax1.add_patch(circ)
		ax2 = plt.subplot2grid((1,4), (0,2),colspan=1)
		ax2.set_title(r'$\textbf{{control}}$')
		ax2.set_xticklabels([])
		ax2.set_yticklabels([])
		ax2.set_adjustable('box-forced')
		dat = result_stack[2][0]
		protein_centers = result_stack[2][1]
		img = ax2.imshow(array(dat).T,interpolation='nearest',origin='LowerLeft',vmax=extremum,vmin=-extremum,
			cmap='bwr',extent=[0,numgridpts,0,numgridpts])
		nprots = 0
		if nprots > 0:
			for pt in protein_centers:
				circ = plt.Circle((int(round(pt[0]/vecs[0]*numgridpts)),
					int(round(pt[1]/vecs[0]*numgridpts))),radius=1.5/64*numgridpts,color='k',
					alpha=1.0)
				ax2.add_patch(circ)
		#divider = make_axes_locatable(ax2)
		#cax = divider.append_axes("right", size="5%", pad=0.05)
		#plt.colorbar(img,ax=ax2,cax=cax)
		#cbar_ax = fig.add_axes([0.85, 0.4, 0.03, 0.75])
		cax = inset_axes(ax2,
		             width="5%",
		             height="100%",
		             bbox_transform=ax2.transAxes,
		             bbox_to_anchor=(0.3, 0.1, 1.05, 0.95),
		             loc= 1)
		fig.colorbar(img,cax=cax)
		cax.tick_params(labelsize=10) 
		cax.set_ylabel(r'$\mathsf{C_{0}(nm^{-1})}$',fontsize=10)
		#fig.subplots_adjust(left=0.05, bottom=0.2, right=0.8, top=0.95, wspace=0.2, hspace=0.2) 
		imagefile = storedir+'/stressmap'+\
			'-span'+str(span)+\
			'-nnum'+str(nnum)+\
			'-ngp'+str(numgridpts)+\
			'-df'+str(distance_factor)+\
			'.png'
		plt.tight_layout() 
		plt.savefig(imagefile,dpi=500,bbox_inches='tight')
		if plotvisible:
			plt.show()
		else:
			plt.clf()
		print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'

#---Run a large parameter sweep of each system individually
if batch_parameter_sweep:
	starttime = time.time()
	print 'Starting analysis job.'
	#---Combine the parameters for sweeping
	ordering = [scan_span,scan_nnum,scan_numgridpts,scan_distance_factor]
	tests = [[i,j,k,l] for i in ordering[0] for j in ordering[1] for k in ordering[2] for l in ordering[3]]
	#---Sweep parameters and render plots
	for ad in analysis_descriptors[analysis_plan]:
		#---Load global variables with calculation specifications used in analysis functions above.
		(systemname,msetfile,picklefile,nprots) = ad
		print 'Running postprocessing: stress tensor map system '+systemname+'.'
		mset = unpickle(pickles+msetfile)
		file3dpp = pickles+picklefile
		#---Make a directory for figures
		if not os.path.isdir(storedir):
			os.mkdir(storedir)
		logmaxminmean = []
		#---Perform global processing steps for this system
		dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
		mset.calculate_average_surface()
		vecs = mean(mset.vecs,axis=0)
		griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
		rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
		xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
			griddims[1]+2)])
		unzipsurfmean = mset.unzipgrid(mset.surf_mean,vecs=vecs)
		interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],
			fill_value=0.0)
		blocksize = (vecs/griddims)
		themeanpos = mean(mset.surf_position)
		if flatref == False:
			avginterpsurf = array([[round((interp(i,j)+themeanpos)/blocksize[2])
				for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
		else:
			avginterpsurf = array([[themeanpos/blocksize[2]
				for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
		if ((themeanpos/blocksize[2] - span < 0.0) or 
			(themeanpos/blocksize[2] + span >  max(dat3dpp[:,2]))):
			print 'Warning: your span exceeds your box dimensions'
		#---Done global steps, now proceed to sweep analysis parameters
		for t in range(len(tests)):
			print 'Running postprocessing: stress tensor map scan '+str(tests[t])+'.'
			span,nnum,numgridpts,distance_factor = tests[t]
			res = calculate_stressmap(span,nnum,numgridpts,distance_factor,plotshow=False,
				plotunsmoothed=False,brokeversion=True,
				imagefile=storedir+'/stressmap.'+systemname+
				'-span'+str(span)+
				'-nnum'+str(nnum)+
				'-ngp'+str(numgridpts)+
				'-df'+str(distance_factor)+
				'.png')
			plot_stressmap(res[0],protein_centers[0],nprots,imagefile=None,plotvisible=False)
		if erase_when_finished:
			del mset
		fp = open(logfile+'.'+systemname,'w')
		for line in logmaxminmean:
			for item in line:
				fp.write(str(item)+'\t')
			fp.write('\n')
		fp.close()
	print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'

#---Run a large parameter sweep of each system individually
if batch_parameter_sweep_framewise:
	starttime = time.time()
	print 'Starting analysis job.'
	(systemname,msetfile,picklefile,nprots) = analysis_descriptors[2]
	mset = unpickle(pickles+msetfile)
	test = [4,32,64,1]
	logmaxminmean = []
	span,nnum,numgridpts,distance_factor = test
	res_collection = []
	for frame in range(len(mset.surf)):
		print 'running frame = '+str(frame)
		file3dpp = pickles+'/'+datdir3dpp+'/'+'md.part0008.fr'+str('%04d'%frame)+'.lp.dat3d'
		dat3dpp = array([[float(i) for i in line.strip().split()] for line in open(file3dpp)])
		griddims = [int(max(dat3dpp[:,i])) for i in range(3)]
		vecs = mean(mset.vecs,axis=0)
		blocksize = (vecs/griddims)
		unzipsurfmean = mset.unzipgrid(mset.surf[frame],vecs=vecs)
		themeanpos = mset.surf_position[frame]
		rawresults = [[[] for j in range(griddims[1]+1)] for i in range(griddims[0]+1)]
		xypts = array([[i,j] for i in linspace(0,vecs[0],griddims[0]+2) for j in linspace(0,vecs[1],
			griddims[1]+2)])
		interp = scipy.interpolate.LinearNDInterpolator(unzipsurfmean[:,0:2],unzipsurfmean[:,2],
			fill_value=0.0)
		avginterpsurf = array([[round((interp(i,j)+themeanpos)/blocksize[2]) 
			for i in linspace(0,vecs[0],griddims[0]+2)] for j in linspace(0,vecs[1],griddims[1]+2)])
		if ((themeanpos/blocksize[2] - span < 0.0) or 
			(themeanpos/blocksize[2] + span >  max(dat3dpp[:,2]))):
			print 'Warning: your span exceeds your box dimensions'
		res = calculate_stressmap(span,nnum,numgridpts,distance_factor,plotunsmoothed=False,
			brokeversion=False)
		res_collection.append(res)
		#plot_stressmap(res[0],res[1],nprots,numgridpts,imagefile=None,plotvisible=True)
	pickle.dump(res_collection,open('/home/rpb/v550rescol','w'))
	'''
	if erase_when_finished:
		del mset
	fp = open(logfile+'.'+systemname,'w')
	for line in logmaxminmean:
		for item in line:
			fp.write(str(item)+'\t')
		fp.write('\n')
	fp.close()
	'''
	print 'Job complete and it took '+str(1./60*(time.time()-starttime))+' minutes.'


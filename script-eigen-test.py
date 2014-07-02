#!/usr/bin/python

#---input files
basedir = '/'
#grofile = '/home/rpb/compbio-alt/membrane-v612/t2-pip2-fix/md.part0005.gro'
#traj = '/home/rpb/compbio/membrane-v612-enthx1-12800/s9-trestles/md.part0003.10000-80000-200.xtc'
grofile = '/home/rpb/worker/covtest/protein.gro'
traj = '/home/rpb/worker/covtest/protein-short.xtc'

#---custom imports
interact = True
from membrainrunner import *
from mayavi import mlab
execfile('locations.py')
execfile('header-cgmd.py')
universe = Universe(basedir+'/'+grofile,basedir+'/'+traj)

#---store the points
pts = []
nframes = len(universe.trajectory)
import MDAnalysis
import MDAnalysis.analysis
ref="ref.gro"
targetdir = './'
gro_writer=MDAnalysis.coordinates.GRO.GROWriter(targetdir+ref)
gro_writer.write(universe,0)
reference=Universe(targetdir+ref).selectAtoms('name CA')

for frameno in range(nframes):
	status('status: storing points for fr = '+str(frameno+1)+'/'+str(nframes))
	frame = universe.trajectory[frameno]
	selection = universe.selectAtoms('name CA')
	MDAnalysis.analysis.align.alignto(selection,reference,mass_weighted=False)
	pts.append(selection.coordinates())



#---compute covariance matrix (a bit slow, but can be optimized)
npts = shape(pts)[1]
avg = mean(pts,axis=0)
avg1 = ravel(avg)
if 'covmat' not in globals():
	covmat = [[[] for j in range(len(avg1))] for i in range(len(avg1))]
	for i in range(len(avg1)):
		for j in range(len(avg1)):
			status('i = '+str(i))
			#covmat[i][j] = mean([(ravel(pts[fr])[i]-avg1[j]) for fr in range(len(pts))])
			covmat[i][j] = mean([(ravel(pts[fr])[i]-avg1[i])*(ravel(pts[fr])[j]-avg1[j]) for fr in range(len(pts))])
		
#---plot covariance matrix
#plt.imshow(array(covmat),interpolation='nearest',cmap=mpl.cm.binary)
#plt.savefig('/home/rpb/covariance.png')
#plt.show()

joecov = pickle.load(open('/home/rpb/worker/covtest/covariance-joe','r'))
rl = []
fp = open('/home/rpb/worker/covtest/protein-short.dcd.varcov.dat','r')
for line in fp:
	rl.append([float(i) for i in line.split()])
rl = array(rl)

fig=plt.figure()
ax=plt.subplot(131)
ax.imshow(array(covmat),interpolation='nearest',cmap=mpl.cm.binary)
ax=plt.subplot(132)
ax.imshow(array(joecov),interpolation='nearest',cmap=mpl.cm.binary)
ax=plt.subplot(133)
ax.imshow(array(rl-covmat),interpolation='nearest',cmap=mpl.cm.binary)
plt.show()



if 0:
	#---compute eigenvalues
	eigenvector_look = 0
	eigvals,eigvecs = numpy.linalg.eig(covmat)
	eigvecs3 = array(abs(eigvecs[:,eigenvector_look])).reshape((npts,3))
	
	#---show average structure (black) alongside one of the modes (white) in 3D
	meshpoints(avg,scale_factor=2)
	meshpoints(avg+50*eigvecs3,scale_factor=2,color=(1,1,1))
	
	#---plot the eigenspectrum
	plt.plot(abs(eigvals))
	plt.yscale('log')
	plt.title('Eigenspectrum')
	plt.show()


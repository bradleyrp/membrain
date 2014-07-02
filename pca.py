from scipy import linalg as la
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.analysis as analysis
import numpy as np
import os
import matplotlib.pylab as plt

class PCA(object):
    """perform principal component analysis"""
    def __init__(self, topology, trajectory, selection='protein and backbone', 
                 targetdir=os.path.curdir):
        self.topology = topology
        self.trajectory = trajectory
        self.targetdir = targetdir
        self.selection = selection
        trajbase = os.path.splitext(os.path.basename(trajectory))[0]
        output = trajbase + '_pca.dat'
        self.output = os.path.join(self.targetdir, output)
        self.universe = Universe(topology, trajectory)
    
    def make_covariance(self):
        selection = self.universe.selectAtoms(self.selection)
        #reference = self.universe.selectAtoms(self.selection)
        ref="ref.gro"
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.targetdir+ref)
        gro_writer.write(self.universe,0)
        reference=Universe(self.targetdir+ref).selectAtoms(self.selection)
        num_atoms=selection.numberOfAtoms()
        num_frames=self.universe.trajectory.numframes
        dof=num_atoms*3
        cov=np.zeros((dof,dof))
        
        print "covariancs matrix will have shape {0} from {1} frames".format(np.shape(cov),num_frames)
        
        #num_confs=0
        coordsum=np.zeros(dof)
        for ts in self.universe.trajectory:
            analysis.align.alignto(selection,reference,mass_weighted=True)
            coords = selection.coordinates().flatten()
            coordsum += coords
            cov += np.outer(coords,coords)

            # mobile = selection.coordinates() - selection.centerOfMass()
            # ref = reference.coordinates() -reference.centerOfMass()
            # R = analysis.align.rotation_matrix(mobile,ref)
            # selection.atoms.translate(-selection.centerOfMass())
            # selection.atoms.rotate(R[0])
            # selection.atoms.translate(reference.centerOfMass())
            # coords = selection.coordinates().flatten()
            # coordsum += coords
            # cov += np.outer(coords,coords)
        import matplotlib
        cov /= num_frames            
        coordsum /= num_frames
        cov -= np.outer(coordsum,coordsum)
        plt.imshow(cov,cmap=matplotlib.cm.binary)
        plt.show()
        return cov


"""    
        masses = np.repeat(selection.masses(), 3)
        mass_matrix = np.sqrt(np.identity(len(masses))*masses)
        #cov1 = np.dot(cov,mass_matrix)
        #self.covariance = np.dot(mass_matrix, cov1)
        cov2=cov.flatten()
        #print cov[0]
        for j in range(1):
            for i in range(0,15,3):
                pass
                #print cov2[dof*j+i], cov2[dof*j+i+1], cov2[dof*j+i+2]
        #eigvals,eigvecs=la.eig(self.covariance)
        #fh = open('cov.dat','w')
        #for var in enumerate(self.covariance):
        #fh.write(self.covariance)
        #fh.close()
        #return self.covariance
"""

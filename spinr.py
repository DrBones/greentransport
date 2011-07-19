# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import scipy
#import matplotlib.pyplot as plt
#import time
#import scipy.linalg as sl
scipy.set_printoptions(precision=3,suppress=True)
from world import World
from model import Model
#from sparseblockslice import SparseBlocks
from aux import spy as sspy
#from io import writeVTK
def main():
    global device, model
    device = World()
    model = Model(device)
    model._Model__build_H()

def alt():
    global smodel, sdevice
    from scipy import asarray
    sdevice = World('canvas/wire70x30spinorbit.bmp')
    smodel = Model(sdevice)
    smodel.block_sizes = asarray([smodel.wafer.shape[1]]*smodel.wafer.shape[0])
    smodel.simpleH()
    smodel.eigensigma()
    print 'Hamiltonian shape = ', smodel.H.shape
    print 'Sigma shape (Mode: Normal) = ', smodel.sigma(smodel.contacts[0], 0.2).shape
    smodel.build_convolutionkernel()

    #smodel.simpleenergyintegrate(smodel.LRGM)

def sweep(instance):
    from scipy import linspace
    from evtk.vtk import VtkGroup
    from io import writeVTK
    g = VtkGroup("./group")
    i=0
    for energy_multi in linspace(-0.001,0.1,500):
        print i
        dens = instance.dorrgm(energy_multi*instance.t0)
        dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        writeVTK('output/huzzah'+str(i), 49, 99, pointData={"Density":dens})
        g.addFile(filepath='output/huzzah'+str(i)+'.vtr', sim_time=i)
        i+=1
    g.save()

def spinint(instance):
    from scipy import linspace,zeros
    #from evtk.vtk import VtkGroup
    from io import writeVTK
    #g = VtkGroup("./spingroup")
    i=0
    dens = zeros((instance.wafer.shape[0],instance.wafer.shape[1]))
    Espace =linspace(instance.potential_drop[1],instance.potential_drop[0],500)
    dE = Espace[1]-Espace[0]
    for energy_multi in Espace:
        print i
        dens = dens + instance.spindens(energy_multi*instance.t0)*dE
        #g.addFile(filepath='output/huzzah'+str(i)+'.vtr', sim_time=i)
        i+=1
    writeVTK('spindens', 29, 69, pointData={"Spin Density":dens.flatten()})
    return dens
    #g.save()

if __name__ == '__main__':
    alt()

#A, sigma_in_l, sigma_in_r = model.build_A(0.1)
#Ablock = SparseBlocks(A, model.block_sizes)
#green_diag, grl, Gr = model.RRGM(Ablock)

#timeittook=time.time()
#def runrrgm():
#    ldos = scipy.array([0]*len(d.nodes))
#    dos = scipy.array([0]*len(d.nodes))
#    for energy in Parameters.Egrid:
#        print energy
#        tic = time.time()
#        sigma_l = c.build_sigma(d.nodes, d.contact[0],energy - Parameters.potential_drop[0])
#        sigma_r = c.build_sigma(d.nodes, d.contact[1], energy - Parameters.potential_drop[1])
#        A = energy*scipy.sparse.eye(len(d.nodes),len(d.nodes),dtype=scipy.complex128, format='lil') - d.HD - sigma_l  - sigma_r
#        block, Ablock = blocks(d.conductor, A)
#        grl, Gr = RRGM(block, Ablock)
#        print time.time() -tic
#        tic = time.time()
#        sigma_in_l = -2* sigma_l.imag[0:block[0],0:block[0]] * fermifunction(energy, mu=Parameters.mu_l)
#        sigma_in_r = -2* sigma_r.imag[-block[-1]:,-block[0]:] * fermifunction(energy, mu=Parameters.mu_r)
#        Gl = LRGM(block, Ablock, grl,Gr, sigma_in_l, sigma_in_r, energy)
#        print time.time() -tic
#        tic = time.time()
#        green_diag = scipy.array([])
#        for i in range(len(block)):
#            diag = scipy.array(Gr[i,i].diagonal()).reshape(-1)
#            green_diag = scipy.hstack((green_diag, diag))
#        dos = dos + fermifunction(energy)*green_diag.imag*Parameters.dE/(-2*scipy.pi*Parameters.a)
#        less_green_diag = scipy.array([])
#        for i in range(len(block)):
#            less_diag = scipy.array(Gl[i,i].diagonal()).reshape(-1)
#            less_green_diag = scipy.hstack((less_green_diag, less_diag))
#        ldos = ldos + less_green_diag.real*Parameters.dE/(scipy.pi*Parameters.a)
#        #dos=dos+(fermifunction(energy)*Parameters.dE.real*green_diag.imag/(-2*scipy.pi*parameters.a))
#        print time.time() -tic
#    return dos,ldos, Gr, Gl
#
#dos,ldos, Gr, Gl = runrrgm()
#
#
##lgmatrix = sl.block_diag(lg[0,0], lg[1,1], lg[2,2], lg[3,3], lg[4,4], lg[5,5], lg[6,6], lg[7,7])
##gmatrix = sl.block_diag(Gr[0,0], Gr[1,1], Gr[2,2], Gr[3,3], Gr[4,4], Gr[5,5], Gr[6,6], Gr[7,7])
##plt.imshow(scipy.absolute(scipy.asarray(d.HD().todense())))
##plt.imshow(scipy.asarray(HD.todense()).real)
##imshow(absolute(lgmatrix))
##imshow(-gmatrix.imag/2*scipy.pi)
##plt.show()
#
#print "The Process took", time.time()-timeittook, "seconds"
##Nodes2 = compose_geo2(Cond)

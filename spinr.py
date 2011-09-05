# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import scipy
#import matplotlib.pyplot as plt
#import time
#import scipy.linalg as sl
#scipy.set_printoptions(precision=3,suppress=True)
from world import World
from model import Model
#from sparseblockslice import SparseBlocks
from aux import spy as sspy
#from io_spinr import writeVTK
def main():
    global device, model
    device = World()
    model = Model(device)
    model._Model__build_H()

def alt():
    global smodel, sdevice
    from scipy import asarray
    sdevice = World('canvas/wire70x30fullspinorbit.bmp')
    smodel = Model(sdevice)
    print 'Bias used: ',smodel.potential_drop
    print 'Fermienergy used: ',smodel.Efermi
    #smodel.block_sizes = asarray([smodel.wafer.shape[1]]*(smodel.wafer.shape[0]))
    #smodel.simpleH()
    #smodel.eigensigma()
    #print 'Hamiltonian shape = ', smodel.H.shape
    #print 'Sigma shape (Mode: Normal) = ', smodel.sigma(smodel.contacts[0], smodel.Efermi).shape
    #smodel.build_convolutionkernel()

    #smodel.simpleenergyintegrate(smodel.LRGM)

def sweep(instance,name=''):
    from scipy import linspace,zeros,array,sum,trace,pi
    from evtk.vtk import VtkGroup
    from io_spinr import writeVTK
    import datetime
    now = datetime.datetime.now()
    name = '{0}.{1}{2}'.format(now.day, now.hour, now.minute)
    from matplotlib.backends.backend_pdf import PdfPages
    #from pylab import plot,figure,title,close,imshow
    import matplotlib.pyplot as plt
    #g = VtkGroup("./group"+name)
    conductivity = []
    i=0
    pdf = PdfPages('Density_and_Conductivity.pdf')
    intdens = zeros((instance.canvas[0],instance.canvas[1]))
    for enmulti in linspace(instance.Efermi-0.09*instance.t0,instance.Efermi+0.5*instance.t0,500):
        print i
        lrgm_val = instance.dolrgm(enmulti)
        filename = 'output/spinr'+name+'_'+str(i)
        if instance.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            writeVTK(filename, 29, 199, pointData={"Density":edens,"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            #writeVTK(filename, 29, 199, pointData={"Density":edens})
        t = instance.transmission(instance.grl)
        G= array(trace(abs(t)**2))
        conductivity.append(G)
        #dens = instance.dorrgm(energy_multi*instance.t0)
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        #g.addFile(filepath=filename+'.vtr', sim_time=i)
        plt.figure(figsize=(3,3))
        plt.imshow(edens)
        plt.colorbar()
        #plot(instance.v.imag)
        #plot(instance.v.real)
        plt.title('Page '+str(i))
        pdf.savefig()
        #close()
        i+=1
    #g.save()
    plt.figure(figsize=(3,3))
    plt.plot(conductivity)
    pdf.savefig()
    pdf.close()
    #import pudb; pudb.set_trace()
    return

def spinint(instance):
    from scipy import linspace,zeros
    #from evtk.vtk import VtkGroup
    from io_spinr import writeVTK
    instance.spinH()
    instance.setmode('spin')
    #g = VtkGroup("./spingroup")
    i=0
    dens = zeros((instance.wafer.shape[0],instance.wafer.shape[1]))
    Espace =linspace(instance.potential_drop[1]/2,instance.potential_drop[0]/2,50)
    dE = Espace[1]-Espace[0]
    for E_rel in Espace:
        print i, 'Rel Energy: ',E_rel,'Total Energy: ',instance.Efermi+E_rel
        dens = dens + instance.spindens(E_rel)*dE
        #g.addFile(filepath='output/huzzah'+str(i)+'.vtr', sim_time=i)
        i+=1
        writeVTK('output/spindens'+str(i), 29, 69, pointData={"Spin Density":dens.imag.flatten()})
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

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
    print "This is spinr, please import as module and supply canvas as: \n spinr.init_with('canvas.bmp')"

def init_with(canvas):
    #global model
    world = World(canvas)
    model = Model(world)
    print 'Size of the world: ', world.canvas
    print 'Fermi energy set: ',model.Efermi
    print 'Bias appied: ', model.potential_drop
    print 'Magnetic field on init: ', model.BField
    return model

def qpc_opening_sweep(instance,name=''):
    from scipy import linspace,zeros,array,sum,trace,pi
    #from evtk.vtk import VtkGroup
    from io_spinr import writeVTK
    import matplotlib
    matplotlib.use('Agg')
    import datetime
    now = datetime.datetime.now()
    name = '{0}.{1}{2}'.format(now.day, now.hour, now.minute)
    from matplotlib.backends.backend_pdf import PdfPages
    #from pylab import plot,figure,title,close,imshow
    import matplotlib.pyplot as plt
    #g = VtkGroup("./group"+name)
    i=0
    pdf = PdfPages('Density_and_Conductivity'+name+'.pdf')
    transmissions = []
    for i in range(200):
        print i
        shift = i/2.5
        print "Setting up Potential Landscape"
        instance.circular_qpc(shift,radius=30,scale=100)
        print "Starting to generate Hamiltonian"
        instance.setmode('normal')
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", '{0}:{1}:{2}'.format(datetime.datetime.now().hour, datetime.datetime.now().minute,datetime.datetime.now().second)
        lrgm_val = instance.dolrgm(instance.Efermi)
        print "Finished lrgm ...Yeah!", '{0}:{1}:{2}'.format(datetime.datetime.now().hour, datetime.datetime.now().minute,datetime.datetime.now().second)
        filename = 'output/spinr'+name+'_'+str(i)
        if instance.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            writeVTK(filename, 29, 199, pointData={"Density":edens,"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            #writeVTK(filename, 29, 199, pointData={"Density":edens})
        t = instance.transmission(instance.grl)
        transmissions.append(t)
        #dens = instance.dorrgm(energy_multi*instance.t0)
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        #g.addFile(filepath=filename+'.vtr', sim_time=i)
        plt.figure(figsize=(3,3))
        plt.imshow(edens)
        plt.colorbar()
        plt.title('Page '+str(i))
        pdf.savefig()
        #close()
        i+=1
    #g.save()
    plt.figure(figsize=(3,3))
    plt.plot(transmissions)
    pdf.savefig()
    pdf.close()
    #import pudb; pudb.set_trace()
    return transmissions

def magnetic_field_sweep(instance,name=''):
    from scipy import linspace
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import datetime
    now = datetime.datetime.now()
    name = '{0}.{1}{2}'.format(now.day, now.hour, now.minute)
    pdf = PdfPages('Density_and_Conductivity'+name+'.pdf')
    transmissions = []
    sampling = 40
    i=0
    for BField in linspace(0,8,sampling):
        print i
        instance.BField = BField
        print 'Magnetic field of',BField,'set'
        print "Starting to generate Hamiltonian"
        instance.setmode('normal')
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", '{0}:{1}:{2}'.format(datetime.datetime.now().hour, datetime.datetime.now().minute,datetime.datetime.now().second)
        lrgm_val = instance.dolrgm(instance.Efermi)
        print "Finished lrgm ...Yeah!", '{0}:{1}:{2}'.format(datetime.datetime.now().hour, datetime.datetime.now().minute,datetime.datetime.now().second)
        #filename = 'output/spinr'+name+'_'+str(i)
        if instance.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            #writeVTK(filename, 29, 199, pointData={"Density":edens,"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            #writeVTK(filename, 29, 199, pointData={"Density":edens})
        t = instance.transmission(instance.grl)
        transmissions.append(t)
        #dens = instance.dorrgm(energy_multi*instance.t0)
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        #g.addFile(filepath=filename+'.vtr', sim_time=i)
        plt.figure(figsize=(4,4))
        plt.imshow(edens)
        plt.title('Electron Density')
        plt.xlabel('Y-Direction')
        plt.ylabel('X-Direction')
        plt.colorbar()
        plt.title('Page '+str(i))
        pdf.savefig()
        #close()
        i+=1
    #g.save()
    plt.figure(figsize=(4,4))
    plt.plot(transmissions)
    plt.title('Transmission')
    plt.xlabel('Magnetic Field')
    plt.ylabel('Transmission')
    pdf.savefig()
    pdf.close()
    #import pudb; pudb.set_trace()
    return transmissions

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
    for enmulti in linspace(instance.Efermi-0.09*instance.t0,instance.Efermi+0.5*instance.t0,70):
        print i
        shift = 100+i
        print "Setting up Potential Landscape"
        instance.circular_qpc(shift,radius=30,scale=100)
        print "Starting to generate Hamiltonian"
        instance.setmode('normal')
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", '{0}:{1}:{2}'.format(now.hour, now.minute,now.second)
        lrgm_val = instance.dolrgm(enmulti)
        print "Finished lrgm ...Yeah!", '{0}:{1}:{2}'.format(now.hour, now.minute,now.second)
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
    main()

# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
#import scipy
#import matplotlib.pyplot as plt
#import time
#import scipy.linalg as sl
#scipy.set_printoptions(precision=3,suppress=True)
from world import World
from model import Model
#from sparseblockslice import SparseBlocks
#from aux import spy as sspy
#from io_spinr import writeVTK
def main():
    print "This is spinr, please import as module and supply canvas as: \n spinr.init_with('canvas.bmp')"

def init_with(canvas=None):
    if canvas is not None:
        world = World(canvas)
        model = Model(world)
        print 'Size of the world: ', model.p.canvas.shape
        print 'Fermi energy set: ',model.p.Efermi
        print 'Bias appied: ', model.p.potential_drop
        print 'Magnetic field on init: ', model.p.BField
        return model
    else:
        print "Please specify atlas as: atlas='filename.bmp'"

def qpc_opening_sweep(instance,name=''):
    from scipy import linspace,zeros,array,sum,trace,pi,complex128
    print 'I am in qpc_opening_sweep now'
    from evtk.vtk import VtkGroup
    from io_spinr import writeVTK
    import gc 
    import numpy as np
    # import matplotlib
    # matplotlib.use('Agg')
    import time
    print 'The time is: ', time.strftime('%X')
    name = time.strftime('%b-%d-%Y-%H.%M')
    import os, errno
    home = os.environ['HOME']
    filepath = home+'/spinr/output/'+str(instance.atlas[-15:-4])+'-'+str(instance.p.task_id)
    try:
        print 'trying to create dir', filepath
        os.makedirs(filepath)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    # from matplotlib.backends.backend_pdf import PdfPages
    #from pylab import plot,figure,title,close,imshow
    # import matplotlib.pyplot as plt
    print 'creating group file'
    g = VtkGroup(filepath+"/group"+name)
    print 'finished creating group file'
    # i=0
    # pdf = PdfPages('Density_and_Conductivity'+name+'.pdf')
    sweep_range = 400
    transmissions = np.zeros(sweep_range,dtype=complex128)
    print 'Setting mode to graph'
    instance.setmode('spin_graph')
    print 'done setting mode to spin_graph'
    # shift = -140
    #charge = 8
    slope = 0
    for i in range(sweep_range):
        print '---------------------------------------------------------'
        print 'Step Number: ',i
        print '---------------------------------------------------------'
        step = 0.05/400
        #charge -=step
        slope +=step
        print "Setting up Potential Landscape"
        instance.p.linearsmooth_qpc(slope,scale=0.56*instance.p.t0,xi=10)
        #instance.p.pointcharge_qpc(charge=charge, scale = 1)
        print "Starting to update Hamiltonian"
        instance.update_hamil_diag()
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", time.strftime('%X')
        lrgm_val = instance.dolrgm(instance.p.Efermi)
        print "Finished lrgm ...Yeah!", time.strftime('%X')
        filename = filepath+'/tstub'+name+'_'+str(i)
        if instance.p.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            writeVTK(filename, instance.p.canvas.shape[1]-1, instance.p.canvas.shape[0]-1, pointData={"Density":edens[0],"UpDensity":edens[1],"DownDensity":edens[2],"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            writeVTK(filename, instance.p.canvas.shape[1]-1,instance.p.canvas.shape[0]-1 , pointData={"Density":edens[0]})
        del lrgm_val
        del edens
        try:
            del spindens
        except NameError:
            pass
        t = instance.transmission(instance.grl)
        transmissions[i]=t
        np.save(filepath+'/'+name+'_transmission',transmissions)
        print 'Calculated transmission: ',t
        del instance.grl
        transmissions[i]=t
        gc.collect()
        #dens = instance.dorrgm(energy_multi*instance.t0)
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        g.addFile(filepath=filename+'.vtr', sim_time=i)
        # plt.figure(figsize=(3,3))
        # plt.imshow(edens)
        # plt.colorbar()
        # plt.title('Page '+str(i))
        # pdf.savefig()
        #close()
        # i+=1
    g.save()
    # plt.figure(figsize=(3,3))
    # plt.plot(transmissions)
    # pdf.savefig()
    # pdf.close()
    #import pudb; pudb.set_trace()
    return transmissions

def energy_sweep(instance,name=''):
    from scipy import linspace,zeros,array,sum,trace,pi,sin,cos,complex128
    print 'I am in energy_sweep now'
    from evtk.vtk import VtkGroup
    import numpy as np
    from io_spinr import writeVTK
    import gc 
    # import matplotlib
    # matplotlib.use('Agg')
    import time
    print 'The time is: ', time.strftime('%X')
    name =str(instance.atlas[-15:-4])
    creation_time = str(time.strftime('%b-%d-%Y-%H.%M'))
    import os, errno
    home = os.environ['HOME']
    filepath = home+'/spinr/output/'+name+'-'+creation_time+'-'+str(instance.p.task_id)+'/'
    try:
        print 'trying to create dir', filepath
        os.makedirs(filepath)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    # from matplotlib.backends.backend_pdf import PdfPages
    #from pylab import plot,figure,title,close,imshow
    # import matplotlib.pyplot as plt
    print 'creating group file'
    g = VtkGroup(filepath+name)
    print 'finished creating group file'
    # i=0
    # pdf = PdfPages('Density_and_Conductivity'+name+'.pdf')
    sweep_range = 100
    transmissions = np.zeros(sweep_range,dtype=complex128)
    print 'Setting mode to spin_graph'
    instance.setmode('spin_graph')
    print 'done setting mode to spin_graph'
    #shift = -140
    #charge = 8
    El=2*instance.p.t0*(1-cos(pi/(instance.p.canvas.shape[1]-2)))
    energy = linspace(6.0/100*40*instance.p.t0/1000,12.0/100*40*instance.p.t0/1000,sweep_range)
    instance.p.potential_drop = [0.004*instance.p.t0/2,-0.004*instance.p.t0/2]
    lines = ['Fermi energy: '+str(instance.p.Efermi)+'\n','Grid spacing: '+str(instance.p.a)+'\n','Hopping parameter t0:'+str(instance.p.t0)+'\n','Spin hopping parameter tSO:'+str(instance.p.tso)+'\n','Potential drop:'+str(instance.p.potential_drop)+'\n','Effective mass:'+str(instance.p.mass)+'\n','Energy range swept: '+str(energy)+'\n','Calculated first mode energy: '+str(El)+'\n']
    with open(filepath+'summary.txt','w') as f:
        f.writelines(lines)
    for i in range(sweep_range):
        print '---------------------------------------------------------'
        print 'Step Number: ',i
        print '---------------------------------------------------------'
        #shift += i
        #step = 8.0/400
        #charge -=step
        #print "Setting up Potential Landscape"
        #instance.p.rectangular_qpc(shift,width=100,scale=100)
        #instance.p.pointcharge_qpc(charge=charge, scale = 1)
        instance.p.stepgrid(20,20)
        print "Starting to update Hamiltonian"
        instance.update_hamil_diag()
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", time.strftime('%X')
        lrgm_val = instance.dolrgm(energy[i])
        print "Finished lrgm ...Yeah!", time.strftime('%X')
        filename = filepath+'/'+name+'_'+str(i)
        if instance.p.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            np.save('filename'+'edens', edens)
            np.save('filename'+'spindens', spindens)
            writeVTK(filename, instance.p.canvas.shape[1]-1, instance.p.canvas.shape[0]-1, pointData={"Density":edens[0],"UpDensity":edens[1],"DownDensity":edens[2],"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            writeVTK(filename, instance.p.canvas.shape[1]-1,instance.p.canvas.shape[0]-1 , pointData={"Density":edens[0]})
        del lrgm_val
        del edens
        try:
            del spindens
        except NameError:
            pass
        t = instance.transmission(instance.grl)
        transmissions[i]=t
        np.save(filepath+'/'+name+'_transmission',transmissions)
        print 'Calculated transmission: ',t
        del instance.grl
        transmissions[i]=t
        gc.collect()
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        g.addFile(filepath=filename+'.vtr', sim_time=i)
        # i+=1
    g.save()
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

def sweep(instance,sweep_range,sweep_type,sweep_min,sweep_max,mode='spin_graph',name=''):
    from scipy import linspace,zeros,array,sum,trace,pi,sin,cos,complex128
    print 'I am in energy_sweep now'
    from evtk.vtk import VtkGroup
    import numpy as np
    from io_spinr import writeVTK
    import gc 
    # import matplotlib
    # matplotlib.use('Agg')
    import time
    print 'The time is: ', time.strftime('%X')
    name =str(instance.atlas[-15:-4])
    if 'creation_time' in dir(instance.p):
        creation_time=instance.p.creation_time
    else:
        creation_time = str(time.strftime('%b-%d-%Y-%H.%M'))
    import os, errno
    home = os.environ['HOME']
    if len(str(instance.p.task_id)) <1:
        filepath = home+'/spinr/output/'+str(instance.p.job_id)+'-'+name+'-'+creation_time+'/'
    else:
        filepath = home+'/spinr/output/'+str(instance.p.job_id)+'-'+name+'-'+creation_time+'/'+str(instance.p.task_id)+'/'
    try:
        print 'trying to create dir', filepath
        os.makedirs(filepath)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    # from matplotlib.backends.backend_pdf import PdfPages
    #from pylab import plot,figure,title,close,imshow
    # import matplotlib.pyplot as plt
    print 'creating group file'
    g = VtkGroup(filepath+name)
    print 'finished creating group file'
    # i=0
    # pdf = PdfPages('Density_and_Conductivity'+name+'.pdf')
    transmissions = np.zeros(sweep_range,dtype=complex128)
    print 'Setting mode to: ',mode
    instance.setmode(mode)
    print 'done setting mode'
    #shift = -140
    charge_range = linspace(sweep_min,sweep_max,sweep_range)
    if sweep_type == 'energy':
        energy = linspace(sweep_min,sweep_max,sweep_range)
    else:
        energy = [instance.p.energy]*sweep_range
        qpc_range=linspace(sweep_min,sweep_max,sweep_range)
    instance.p.energy = energy
    savestats(instance,mode,sweep_type,filepath)
    # instance.p.potential_drop = [0.004*instance.p.t0/2,-0.004*instance.p.t0/2]
    #slope_range=linspace(0.24,0,sweep_range)
    for i in range(sweep_range):
        print '---------------------------------------------------------'
        print 'Step Number: ',i
        print '---------------------------------------------------------'
        #shift += i
        # charge -=step
        print "Setting up Potential Landscape"
        if sweep_type == 'qpccircular':
            instance.p.circular_qpc(shift=qpc_range[i], radius =30, scale=100)
        elif sweep_type == 'qpcpoint':
            instance.p.pointcharge_qpc(charge=charge_range[i], scale = 1)
        elif sweep_type == 'qpcvariational':
            instance.p.linearsmooth_qpc(width=qpc_range[i],slope=instance.p.slope,scale=0.56*instance.p.t0,xi=10)
        elif sweep_type == 'qpcrect':
            instance.p.rectangular_qpc(shift=qpc_range[i],scale=100, width=30)
        elif sweep_type == 'qpctriangular':
            instance.p.triangular_qpc(shift=qpc_range[i],width=50,radius=100,scale=100)
        elif sweep_type == 'ring':
            pass
            #shift = 0 is closed, shift=200 is open
        #elif sweep_type == 'energy':
        #    instance.p.stepgrid(4,4)
        print "Starting to update Hamiltonian"
        instance.update_hamil_diag()
        print "Hamiltonian set up, calculating lrgm (crunch...crunch)", time.strftime('%X')
        lrgm_val = instance.dolrgm(energy[i])
        print "Finished lrgm ...Yeah!", time.strftime('%X')
        filename = filepath+'/'+name+'_'+str(i)
        if instance.p.multi == 2:
            edens =instance.edens(lrgm_val)
            spindens =instance.spindens(lrgm_val)
            np.save(filename+'edens', edens)
            np.save(filename+'spindens', spindens)
            writeVTK(filename, instance.p.canvas.shape[1]-1, instance.p.canvas.shape[0]-1, pointData={"Density":edens[0],"UpDensity":edens[1],"DownDensity":edens[2],"SpinDensity":spindens})
        else:
            edens =instance.edens(lrgm_val)
            np.save(filename+'edens', edens)
            writeVTK(filename, instance.p.canvas.shape[1]-1,instance.p.canvas.shape[0]-1 , pointData={"Density":edens[0]})
        del lrgm_val
        del edens
        try:
            del spindens
        except NameError:
            pass
        t = instance.transmission(instance.grl)
        transmissions[i]=t
        np.save(filepath+'/'+name+'_transmission',transmissions)
        print 'Calculated transmission: ',t
        del instance.grl
        transmissions[i]=t
        gc.collect()
        #dens = -dens.imag/(instance.a**2)*instance.fermifunction(energy_multi*instance.t0, instance.mu)
        #intdens = intdens + spindens
        g.addFile(filepath=filename+'.vtr', sim_time=i)
        # i+=1
    g.save()
    return transmissions

def savestats(instance,mode,sweep_type,filepath):
    lines = ['Fermi energy: '+str(instance.p.energy/instance.p.Efermi)+'\n','Grid spacing: '+str(instance.p.a)+'\n','Hopping parameter t0:'+str(instance.p.t0)+'\n','Spin hopping parameter tSO:'+str(instance.p.tso)+'\n','Potential drop:'+str(instance.p.potential_drop)+'\n','Effective mass:'+str(instance.p.mass)+'\n','Energy range swept: '+str(instance.p.energy)+'\n','Calculated first mode energy: '+str(instance.p.El)+'\n','Graph mode: '+mode+'\n','Sweep type: '+sweep_type+'\n']
    with open(filepath+'summary.txt','w') as f:
        f.writelines(lines)
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

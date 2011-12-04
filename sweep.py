print 'Now comes the sweep script'
print 'importing os'
import os
# print 'importing numpy'
import numpy as n
import glob as g
print 'importing spinr'
import spinr
print 'importing time'
import time
home =os.environ['HOME']
# sm = spinr.init_with(home+'/spinr/canvas/wire400x200.bmp')
sm = spinr.init_with(home+'/spinr/canvas/200x200ring_neg90degrees.bmp')
# sm = spinr.init_with(home+'/spinr/canvas/wire200x100.bmp')
# sm = spinr.init_with(home+'/spinr/canvas/wire300x80.bmp')
if 'SGE_TASK_ID' in os.environ:
    task_id_str = os.environ['SGE_TASK_ID']
    try:
        task_id = int(task_id_str)-1
    except ValueError:
        task_id = ''
else:
    task_id = ''
job_id_str = os.environ['SGE_JOB_ID']
job_id = int(job_id_str)
#--------ring sweeeeep --------------------
# ring_list = g.glob(home+'/spinr/canvas/*ring*')
# ring_list.sort(lambda a,b:cmp(int(a.split('_')[2][:-11]),int(b.split('_')[2][:-11])))
# sm = spinr.init_with(ring_list[task_id])
#--------ring sweeeeep  end--------------------
sm.p.task_id = task_id
sm.p.job_id = job_id
if 'JOB_CREATION_TIME' in os.environ:
    print 'using one time'
    sm.p.creation_time = os.environ['JOB_CREATION_TIME']
print 'This is task: ',task_id,'script staring at: ',time.strftime('%X'),'in',home
# parameter_space = n.linspace(0,1,50)
slope_range=n.linspace(0,0.4,50)
# sm = spinr.init_with(home+'/spinr/canvas/wire200x100.bmp')
# sm.p.Efermi = parameter_space[task_id]*sm.p.Efermi
sm.p.energy=0.16*sm.p.Efermi
print 'Fermi Energy is: ',sm.p.Efermi,'eV (i believe)'
# transmission = spinr.qpc_opening_sweep(sm)
# transmission = spinr.energy_sweep(sm)
# transmission = spinr.sweep(sm,100,'energy',sm.p.El-0.1*sm.p.El,sm.p.El+0.1*sm.p.El,'graph')
# sm.p.linearsmooth_qpc(slope_range[task_id],scale=0.56*sm.p.t0,xi=10)
# sm.p.slope = slope_range[task_id]
# transmission = spinr.sweep(sm,200,'energy',0,0.02*sm.p.Efermi,'spin_graph')
# transmission = spinr.sweep(sm,200,'qpcrect',0,200,'spin_graph')
#transmission = spinr.sweep(sm,200,'qpctriangular',0,200,'graph')
# n.save(home+'/spinr/output/tstub-'+str(task_id)+'/transmission_tstub',transmission)
transmission = spinr.sweep(sm,200,'energy',0,sm.p.Efermi,mode='spin_graph')
# transmission = spinr.sweep(sm,100,'qpcvariational',0,200,mode='graph')
# transmission = spinr.sweep(sm,200,'qpcpoint',40,0,mode='graph')
# transmission = spinr.sweep(sm,200,'qpccircular',72,-30,mode='spin_graph')

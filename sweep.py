print 'Now comes the sweep script'
print 'importing os'
import os
# print 'importing numpy'
# import numpy as n
print 'importing spinr'
import spinr
print 'importing time'
import time
home =os.environ['HOME']
sm = spinr.init_with(home+'/spinr/canvas/wire300x80.bmp')
if 'SGE_TASK_ID' in os.environ:
    task_id_str = os.environ['SGE_TASK_ID']
    try:
        task_id = int(task_id_str)-1
    except ValueError:
        task_id = ''
else:
    task_id = ''
sm.p.task_id = task_id
print 'This is task: ',task_id,'script staring at: ',time.strftime('%X'),'in',home
#parameter_space = n.linspace(0.1,0.2,50)
# sm = spinr.init_with(home+'/spinr/canvas/wire200x100.bmp')
#sm.p.Efermi = parameter_space[task_id]*sm.p.Efermi
# print 'Fermi Energy is: ',sm.p.Efermi,'eV (i believe)'
# transmission = spinr.qpc_opening_sweep(sm)
# transmission = spinr.energy_sweep(sm)
# transmission = spinr.sweep(sm,100,'energy',sm.p.El-0.1*sm.p.El,sm.p.El+0.1*sm.p.El,'graph')
transmission = spinr.sweep(sm,100,'qpc',sm.p.El-0.1*sm.p.El,sm.p.El+0.1*sm.p.El,'graph')
# n.save(home+'/spinr/output/tstub-'+str(task_id)+'/transmission_tstub',transmission)

print 'Now comes the sweep script'
print 'importing os'
import os
print 'importing numpy'
import numpy as n
print 'importing spinr'
import spinr
print 'importing time'
import time
home =os.environ['HOME']
task_id = int(os.environ['SGE_TASK_ID'])-1
print 'This is task: ',task_id,'script staring at: ',time.strftime('%X'),'in',home
parameter_space = n.linspace(0.1,0.2,100)
sm = spinr.init_with(home+'/spinr/canvas/tstub100ax100a.bmp')
sm.p.task_id = task_id
sm.p.Efermi = parameter_space[task_id]*sm.p.Efermi
print 'Fermi Energy is: ',sm.p.Efermi,'eV (i believe)'
transmission = spinr.qpc_opening_sweep(sm)
n.save(home+'/spinr/output/tstub-'+str(task_id)+'/transmission_tstub',transmission)

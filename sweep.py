import numpy as n
import spinr
sm = spinr.init_with('canvas/tstub100ax100a.bmp')
transmission = spinr.qpc_opening_sweep(sm)
n.save('transmission_tstub',transmission)

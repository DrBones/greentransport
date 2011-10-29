import numpy as n
import spinr
sm = spinr.init_with('canvas/200x400wire_template.bmp')
transmission = spinr.qpc_opening_sweep(sm)
n.save('transmission_pointcharge',transmission)

import numpy as np
import matplotlib.pylab as pl
x = np.r_[1:101]

qpcwees = np.load('qpcconductancewees.npy')
qpcme = np.load('qpcme.npy')
#fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'pdf',
          'font.family' : 'serif',
          'font.serif'  :   'computer modern roman',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          # 'text.dvipnghack' : True,
          'text.usetex': True,
          'figure.figsize': fig_size}
pl.rcParams.update(params)
xshift=-0.12
# Plot data
# -------------------------------------
pl.figure()
pl.clf()
pl.axes([0.14,0.2,0.95-0.125,0.95-0.2])
pl.plot(qpcme[:135],label=r'Conductance')
pl.xlabel(r'Width [nm]')
pl.ylabel(r'Conductance [$2e^2/h$]')
pl.yticks(np.linspace(10,0,11),('10','9','8','7','6','5','4','3','2','1','0'))
pl.ylim([0,10])
pl.grid()
#ax = pl.gca()
#ax.yaxis.set_label_coords(xshift, 0.5)
#pl.legend()
pl.savefig('images/qpcme.pdf',transparent='true')
# -------------------------------------
fig = pl.figure()
#ax = fig.add_subplot(111)
#pl.clf()
#fig.clf()
ax=pl.axes([0.14,0.2,0.95-0.125,0.95-0.2])
cax = ax.imshow(qpcwees[10:670,:],interpolation='bicubic',aspect='auto')
#locs, labels = ax.yticks()
pl.yticks(np.linspace(0,600,10),('10','9','8','7','6','5','4','3','2','1','0'))
#pl.yticks(np.linspace(0,600,5),('10','8','6','4','2'))
pl.xticks(np.linspace(0,1100,7),('-2.1','-2.0','-1.8','-1.6','-1.4','-1.2','-1.0'))
#ax = pl.gca()
#ax.xaxis.set_label_coords(0.5, -0.05)
pl.xlabel(r'Gate voltage [V]')
pl.ylabel(r'Conductance [$2e^2/h$]')
pl.grid()
#pl.text(0.85, -0.55, r'$n\times 10^{17}$', fontsize=8, transform = pl.gca().transAxes
#pl.colorbar(orientation='horizontal')
#pl.colorbar()
#pl.text(xshift, 0.85, r'$n$[a.u.]', fontsize=10, transform = pl.gca().transAxes)
fig.savefig('images/qpcwees.pdf',transparent='true')

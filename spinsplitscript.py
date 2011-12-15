import spinr
import numpy as np
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
x = np.r_[1:101]

spindens = np.load('spinsplitkato.npy')
spindenskato = np.load('rawkato.npy')
fig_width_pt = 448.13095  # Get this from LaTeX using \showthe\columnwidth
fig_height_pt = 130  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_height_pt*inches_per_pt  # width in inches
#fig_height = fig_width*golden_mean      # height in inches
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
xshift=1.2
# Plot data
# -------------------------------------
fig = pl.figure()
#ax = fig.add_subplot(111)
#pl.clf()
ax = pl.axes([0.07,0.15,0.95-0.1,0.95-0.1])
im = ax.imshow(spindens,cmap='RdBu_r')
pl.yticks([0,20,40,60,80])
pl.xlabel('$y$\,[$\mu$m]')
pl.ylabel('$x$\,[$\mu$m]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
#pl.yticks([0,20,40,60,80])
#ax = pl.gca()
#ax.xaxis.set_label_coords(0.5, -0.05)
#pl.text(0.85, -0.55, r'$n\times 10^{17}$', fontsize=8, transform = pl.gca().transAxes)
#pl.colorbar(orientation='horizontal')
#pl.colorbar()
cbar = pl.colorbar(im, cax,ticks=[-0.3,0,0.3])
cbar.ax.set_yticklabels(['-1','0','1'])
pl.text(xshift, 0.85, r'$n$ [a.u.]', fontsize=10, transform = pl.gca().transAxes)
fig.savefig('images/spinsplitme.pdf',transparent='true')
# -------------------------------------
fig_width_pt = 448.13095  # Get this from LaTeX using \showthe\columnwidth
fig_height_pt = 150  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_height_pt*inches_per_pt  # width in inches
#fig_height = fig_width*golden_mean      # height in inches
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
fig = pl.figure()
#ax = fig.add_subplot(111)
#pl.clf()
#ax = pl.gca()
#ax.xaxis.set_label_coords(0.5, -0.05)
ax = pl.axes([0.07,0.15,0.95-0.1,0.95-0.1])
im = ax.imshow(spindens,cmap='RdBu_r')
pl.xlabel('$y$\,[$\mu$m]')
pl.ylabel('$x$\,[$\mu$m]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
pl.text(xshift, 0.85, r'$n$ [a.u.]', fontsize=10, transform = pl.gca().transAxes)
#pl.colorbar(orientation='horizontal')
#pl.colorbar()
cbar = pl.colorbar(im, cax,ticks=[-0.3,0,0.3])
cbar.ax.set_yticklabels(['-1','0','1'])
ax.imshow(spindenskato)
fig.savefig('images/spinsplitkato.pdf',transparent='true')

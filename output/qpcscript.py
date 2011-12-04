import numpy as np
import matplotlib.pylab as pl
import glob as g
from mpl_toolkits.axes_grid1 import make_axes_locatable
x = np.r_[1:101]
xshift=1.4
pathlist = g.glob('*wire*')
# Plot data
filepath = '/afs/physnet.uni-hamburg.de/users/ap_n/jsiegl/spinr/images/'
def plotdens(path):
# fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
    fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
    fig_height_pt = 280  # Get this from LaTeX using \showthe\columnwidth
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
    gedens = g.glob(path+'/*edens.npy')
    if len(gedens)<1:
        print 'No edens found'
    else:
        gedens.sort(lambda a,b:cmp(int(a.split('_')[1][:-9]),int(b.split('_')[1][:-9])))
        edens = np.load(gedens[len(gedens)/2])
        pl.figure()
        ax = pl.axes([0,0.1,0.95-0.05,0.95-0.1])
        pl.imshow(edens[0])
        # im = ax.imshow(edens[0])
        # divider = make_axes_locatable(ax)
        pl.xlabel('$y$ [nm]')
        pl.ylabel('$x$ [nm]')
        pl.text(xshift, 0.5, r'electron density $n$ [1/$m^2$]', fontsize=10,verticalalignment='center',rotation='vertical', transform = pl.gca().transAxes)
        # cax = divider.append_axes("top", size="5%", pad=0.05)
        # tickrange = np.linspace(edens.min(),edens.max(),5)
        # label = map(lambda x: r"%.1f" % x, tickrange/1e17)
        pl.colorbar()
        # cbar = pl.colorbar(im, cax=cax,orientation='horizontal')
        # cbar = pl.colorbar(im, cax=cax,orientation='horizontal',ticks=tickrange)
        # cbar.ax.xaxis.set_ticks_position('top')
        # cbar.ax.set_xticklabels(label)
        pl.savefig(filepath+path+'edens.pdf',transparent='true')
# -------------------------------------
def plotspindens(path):
    pl.clf()
    gspindens = g.glob(path+'/*spindens.npy')
    if len(gspindens)<1:
        print 'No spindens.npy found'
    else:
        gspindens.sort(lambda a,b:cmp(int(a.split('_')[1][:-12]),int(b.split('_')[1][:-12])))
        spindens = np.load(gspindens[len(gspindens)/2])
# spindens = spindens[0]
        ax = pl.axes([0,0.1,0.95-0.05,0.95-0.1])
        pl.imshow(spindens,cmap='RdBu_r')
        # divider = make_axes_locatable(ax)
        pl.xlabel('$y$ [nm]')
        pl.ylabel('$x$ [nm]')
        pl.text(xshift, 0.5, r'spin density $s_z$ [$\hbar/m^2$]', fontsize=10,verticalalignment='center',rotation='vertical', transform = pl.gca().transAxes)
        # pl.text(0.5, 1.15, r'$s_z$', fontsize=10, horizontalalignment='center', transform = pl.gca().transAxes)
        pl.colorbar()
        # cax = divider.append_axes("top", size="5%", pad=0.05)
        # tickrange = np.linspace(spindens.min(),spindens.max(),5)
        # label = map(lambda x: r"%.1f" % x, tickrange)
        # cbar = pl.colorbar(im, cax=cax,orientation='horizontal',ticks=tickrange)
        # cbar.ax.xaxis.set_ticks_position('top')
        # cbar.ax.set_xticklabels(label)
        pl.savefig(filepath+path+'spindens.pdf',transparent='true')
# -------------------------------------
def plottransmission(path):
    from matplotlib.ticker import MultipleLocator
# fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
    fig_width_pt = 448.13095/2  # Get this from LaTeX using \showthe\columnwidth
# fig_height_pt = 350  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
# fig_height = fig_height_pt*inches_per_pt  # width in inches
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
    pl.figure()
    pl.clf()
    gtrans = g.glob(path+'/*transmission.npy')
    if len(gtrans)<1:
        print 'No transmissions.npy found'
    else:
        trans = np.load(gtrans[0])
        ax = pl.axes([0.17,0.2,0.95-0.2,0.95-0.2])
        pl.plot(trans)
        pl.xlabel('Shift [points]')
        pl.ylabel('Conductance [$e^2/h$]')
        majorLocator = MultipleLocator(2)
        ax.yaxis.set_major_locator(majorLocator)
        pl.grid()
        pl.savefig(filepath+path+'trans.pdf',transparent=True)

for path in pathlist:
    print path
    plotdens(path)
    plotspindens(path)
    plottransmission(path)

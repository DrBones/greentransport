import spinr
import numpy as np
import matplotlib.pylab as pl
x = np.r_[1:101]

def psi(x,n):
    return np.sqrt(2.0/(len(x)+1))*np.sin(n*np.pi*x/(len(x)+1))
def El(n):
    return 2*sm.p.t0*(1-np.cos((n*np.pi)/100))

# sm = spinr.init_with('canvas/wire200x100.bmp')
# sm.setmode('graph')
# dens2 = sm.edens(sm.dolrgm(El(2)-0.019*El(2)))[0]# maximun near El-0.019*El
# np.save('dens21',dens1)
dens2 = np.load('dens2.npy')
dens1=np.load('dens21.npy')
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
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          # 'text.dvipnghack' : True,
          'text.usetex': True,
          'figure.figsize': fig_size}
pl.rcParams.update(params)
y1 = psi(x,2)/(psi(x,2).max())
y2 = y1**2
y3 = dens2[100,x]
y5 = dens1[100,x]
y4 = (y2*y3.max()+y5)/(y2*y3.max()+y5).max()*y3.max()
xshift=-0.12
# Plot data
# -------------------------------------
pl.figure()
pl.clf()
pl.axes([0.15,0.2,0.95-0.15,0.95-0.2])
pl.plot(x,y1,'g:',label=r'$\Psi$')
pl.plot(x,y2,'-b',label=r'$\Psi^+\Psi$')
pl.xlabel(r'$y$\,[nm]')
pl.ylabel(r'$\pi$\,[arb. units]')
ax = pl.gca()
ax.yaxis.set_label_coords(xshift, 0.5)
pl.legend()
pl.savefig('images/analytical2.pdf',transparent='true')
# -------------------------------------
pl.figure()
pl.clf()
pl.axes([0.15,0.2,0.95-0.15,0.95-0.2])
pl.plot(x,y4,'-b',label='$n2_{analytical}$')
pl.plot(x,y3,'g2',label='$n2_{simulated}$',markersize=3)
pl.plot(x,y5,'r2',label='$n1_{simulated}$',markersize=3)
pl.xlabel(r'$y$\,[nm]')
pl.ylabel(r'$n$\,[1/m$^2$]')
ax = pl.gca()
ax.yaxis.set_label_coords(xshift, 0.5)
locs,labels = pl.yticks()
pl.yticks(locs, map(lambda x: r"$ %.1f $" % x, locs*1e-17))
pl.text(-0.17, 0.9, r'$\times 10^{17}$', fontsize=8, transform = pl.gca().transAxes)
pl.legend()
pl.savefig('images/overlay2.pdf',transparent='true')
# -------------------------------------
pl.figure()
pl.clf()
pl.axes([0.15,0.1,0.95-0.15,0.95-0.1])
pl.imshow(dens2/1e17,aspect='auto')
ax = pl.gca()
ax.xaxis.set_label_coords(0.5, -0.05)
pl.xlabel('$y$\,[nm]')
pl.ylabel('$x$\,[nm]')
pl.text(0.85, -0.55, r'$n\times 10^{17}$', fontsize=8, transform = pl.gca().transAxes)
pl.colorbar(orientation='horizontal')
pl.savefig('images/dens2.pdf',transparent='true')
# -------------------------------------
pl.figure()
pl.clf()
error=(y4-y3)*100/(sum(y3)/len(x))
pl.axes([0.15,0.2,0.95-0.15,0.95-0.2])
pl.plot(x,error,'-r',label=r'Error')
pl.xlabel(r'$y$\,[nm]')
pl.ylabel(r'relative error\,[\%]')
ax = pl.gca()
ax.yaxis.set_label_coords(xshift, 0.5)
locs,labels = pl.yticks()
pl.yticks(locs, map(lambda x: r"$ %.1f $" % x, locs*1e1))
pl.text(-0.19, 0.9, r'$\times 10^{-1}$', fontsize=8, transform = pl.gca().transAxes)
pl.legend()
pl.savefig('images/error2.pdf',transparent='true')

import numpy as np
import pynbody.plot.sph as sph
import pynbody
import matplotlib.pylab as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as clrs
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['figure.figsize'  ] = (3.3,2.55)
plt.rcParams['font.size'       ] = 8
plt.rcParams['legend.fontsize' ] = 8
plt.rcParams['legend.frameon'  ] = False
plt.rcParams['font.family'     ] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction' ] = 'in'
plt.rcParams['ytick.direction' ] = 'in'
plt.rcParams['xtick.top'       ] = True
plt.rcParams['ytick.right'     ] = True
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.75
plt.rcParams['ytick.minor.width'] = 0.5

width = 1500 # Unit: kpc

df = pynbody.load('040/snapshot_030.hdf5')

dens = sph.image(df.g,qty="rho",units="g cm^-3",width=width,cmap="viridis", av_z='rho')

time_snapshot = df.properties['time'].in_units('Gyr') # Unit: Gyr

fig, ax = plt.subplots(nrows=1, ncols=1)

# Here we can get the density normalization
density_norm     =  clrs.LogNorm(5e-26, 1e-28) # g cm^-3

cmapDens = plt.get_cmap('viridis')
extent = (-width/2000, width/2000, -width/2000, width/2000)


fig, ax = plt.subplots(nrows=1, ncols=1)

density = ax.imshow(dens, cmap=cmapDens, norm=density_norm, extent=extent, zorder=1)

ax.annotate(f'Time = {np.around(time_snapshot, 2)} Gyr', xy=(-0.65,-0.65), color='white', zorder=4, fontsize=8)

density_colorbar_ax = fig.add_axes([0.806, 0.125, 0.025, 0.868])
density_colorbar  = plt.colorbar(density, cax=density_colorbar_ax, pad=0)
density_colorbar.set_label(r'$\log \, \rho \;$ (g cm$^{-3}) $', labelpad=8)

fig.subplots_adjust(left=0.13, bottom=0.125, top=0.99, right=0.81, hspace=0.00, wspace=0.0)

ax.set_xlabel(r'$x$ (kpc)')
ax.set_ylabel(r'$y$ (kpc)')
ax.set_aspect('equal')

plt.savefig('density.pdf')
plt.savefig('density.png', dpi=200)

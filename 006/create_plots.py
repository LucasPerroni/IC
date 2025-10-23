import h5py
import pynbody
import numpy as np
import pynbody.plot.sph as sph
import matplotlib.pylab as plt
import matplotlib.colors as clrs
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

SNAPSHOT_CODE = "0005_0006_2000/"
SNAPSHOT_PATH = "/mnt/d/UFPR/IC/snapshots/" + SNAPSHOT_CODE
IMAGE_PATH = "/home/lucasbondep/ic_astronomia/main/006/plots/" + SNAPSHOT_CODE

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

snapshots = []
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    file = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    snapshots.append(file)

count = 0
# s_init, s_end = 90, 100 # 0005_0006_100
# s_init, s_end = 60, 70 # 0005_0006_1000
s_init, s_end = 40, 50 # 0005_0006_2000
for file in snapshots:
    if (count > s_end) | (count < s_init):
    # if (count < 90):
        count += 1
        continue

    print(f"{bcolors.OKGREEN}Generating frame {count} of {len(lines) - 1:03d}...{bcolors.ENDC}")

    # ----------------------------------------------------------------------------------------------

    plt.rcParams['figure.figsize'  ] = (7, 6)
    plt.rcParams['font.size'       ] = 16
    plt.rcParams['legend.fontsize' ] = 16
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

    width = 1800 # Unit: kpc

    # DATA -----------------------------------------------------------------------------------------
    df = pynbody.load(file)
    df.physical_units()
    
    s = h5py.File(file, 'r')

    time_snapshot = s['Header'].attrs[u'Time'] # Gyr
    # time_snapshot = df.properties['time'].in_units('Gyr') # Unit: Gyr

    mi = 0.6 # Average molecular weight
    Mh = 1.67262192 * 10**(-27) # Proton mass in kg
    df.gas["kT"] = (df.gas["u"] * (2 * mi * Mh) / 3) * 6.241506 * 10**15 * 10**(6)

    # DENSITY --------------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1)

    plt.sca(ax)  # Define o eixo atual para o sph.image
    im_array = sph.image(
        df.gas,
        qty="rho",
        units="g cm^-3",
        width=width,
        cmap="twilight",
        av_z='rho',
        vmin=2e-28, 
        vmax=3e-25,
        show_cbar=False,
        noplot=True
    )

    # Constrói o plot manualmente
    extent = (-width/2, width/2, -width/2, width/2)
    norm = clrs.LogNorm(vmin=1e-28, vmax=5e-26)
    im = ax.imshow(im_array, cmap="twilight", extent=extent, norm=norm)
    ax.invert_xaxis()

    fig.subplots_adjust(left=0.15, bottom=0.075, top=0.925, right=0.81, hspace=0.00, wspace=0.0)
    ax.set_xlabel(r'$x$ (kpc)')
    ax.set_ylabel(r'$y$ (kpc)')
    ax.set_aspect('equal')
    ax.annotate(f'{np.around(time_snapshot, 3)} Gyr', xy=(-width/3.5, width/2.3), 
                color='black', zorder=4, fontsize=16)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.1)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(r'$\log \, \rho \;$ (g cm$^{-3}) $', labelpad=8)

    plt.savefig(f'{IMAGE_PATH}density_{count:03d}.png', dpi=300)

    # TEMPERATURE ----------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1)

    plt.sca(ax)  # Define o eixo atual para o sph.image
    im_array = sph.image(
        df.gas,
        qty="kT",
        width=width,
        cmap="inferno",
        show_cbar=False,
        noplot=True
    )

    # Constrói o plot manualmente
    extent = (-width/2, width/2, -width/2, width/2)
    norm = clrs.Normalize(vmin=0, vmax=12)
    im = ax.imshow(im_array, cmap="inferno", extent=extent, norm=norm)
    ax.invert_xaxis()

    fig.subplots_adjust(left=0.15, bottom=0.075, top=0.925, right=0.81, hspace=0.00, wspace=0.0)
    ax.set_xlabel(r'$x$ (kpc)')
    ax.set_ylabel(r'$y$ (kpc)')
    ax.set_aspect('equal')
    ax.annotate(f'{np.around(time_snapshot, 3)} Gyr', xy=(-width/3.5, width/2.3), 
                color='white', zorder=4, fontsize=16)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.1)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(r'$kT$ (keV)', labelpad=8)

    plt.savefig(f'{IMAGE_PATH}temperature_{count:03d}.png', dpi=300)

    count += 1

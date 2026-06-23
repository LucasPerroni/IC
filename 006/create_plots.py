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

LEFT=0.17
BOTTOM=0.055
TOP=0.975
RIGHT=0.81
HSPACE=0.00
WSPACE=0.0
YEAR_X_FRAC=3.5 # 3.5 if DECIMAL_COUNT=1 | 4.5 if DECIMAL_COUNT=3
YEAR_Y_FRAC=2.3
DECIMAL_COUNT=1
DECIMAL_FONT_SIZE=16
SHOW_CBAR=False

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

SNAPSHOT_CONFIG = {
    "0005_0005_0/": {
        "s_init": 80,
        "s_end": 90,
    },
    "0005_0005_1000/": {
        "s_init": 54,
        "s_end": 64,
    },
    "0005_0005_3000/": {
        "s_init": 30,
        "s_end": 40,
    },
    "0005_0006_0/": {
        "s_init": 95,
        # "s_end": 100,
        "s_end": 105,
        "mach_med": 2.18,
    },
    "0005_0006_100/": {
        "s_init": 90,
        # "s_end": 95,
        "s_end": 100,
        "mach_med": 2.24,
    },
    "0005_0006_1000/": {
        "s_init": 60,
        # "s_end": 65,
        "s_end": 70,
        "mach_med": 2.70,
    },
    "0005_0006_2000/": {
        "s_init": 40,
        # "s_end": 45,
        "s_end": 50,
        "mach_med": 2.76,
    },
    "0005_0006_3000/": {
        "s_init": 26,
        # "s_end": 31,
        "s_end": 36,
        "mach_med": 2.62,
    },
    "0005_0007_0/": {
        "s_init": 82,
        "s_end": 92,
    },
    "0005_0007_1000/": {
        "s_init": 53,
        "s_end": 63,
    },
    "0005_0007_3000/": {
        "s_init": 30,
        "s_end": 40,
    },
}

try:
    cfg = SNAPSHOT_CONFIG[SNAPSHOT_CODE]
    s_init = cfg["s_init"]
    s_end = cfg["s_end"]
    if cfg["mach_med"]:
        mach_med = cfg["mach_med"]
except KeyError:
    raise ValueError(f"Configuração não encontrada para SNAPSHOT_CODE = {SNAPSHOT_CODE}")


snapshots = []
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    file = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    snapshots.append(file)

count = 0
for file in snapshots:
    if (count > s_end) | (count < s_init):
        count += 1
        continue

    print(f"{bcolors.OKGREEN}Generating frame {count} of {len(lines) - 1:03d}...{bcolors.ENDC}")

    # ----------------------------------------------------------------------------------------

    plt.rcParams['figure.figsize'  ] = (7, 6)
    plt.rcParams['font.size'       ] = 20
    plt.rcParams['legend.fontsize' ] = 20
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
    # width = 6000 # Unit: kpc

    # DATA ------------------------------------------------------------------------------------
    df = pynbody.load(file)
    df.physical_units()
    
    s = h5py.File(file, 'r')

    time_snapshot = s['Header'].attrs[u'Time'] # Gyr
    # time_snapshot = df.properties['time'].in_units('Gyr') # Unit: Gyr

    mi = 0.6 # Average molecular weight
    Mh = 1.67262192 * 10**(-27) # Proton mass in kg
    df.gas["kT"] = (df.gas["u"] * (2 * mi * Mh) / 3) * 6.241506 * 10**15 * 10**(6)

    # GAS DENSITY -----------------------------------------------------------------------------
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
        show_cbar=SHOW_CBAR,
        noplot=True
    )

    # Constrói o plot manualmente
    extent = (-width/2, width/2, -width/2, width/2)
    norm = clrs.LogNorm(vmin=1e-28, vmax=5e-26)
    im = ax.imshow(im_array, cmap="twilight", extent=extent, norm=norm)
    ax.invert_xaxis()

    fig.subplots_adjust(left=LEFT, bottom=BOTTOM, top=TOP, right=RIGHT, hspace=HSPACE, wspace=WSPACE)
    ax.set_xlabel(r'$x$ (kpc)')
    ax.set_ylabel(r'$y$ (kpc)')
    ax.set_aspect('equal')
    ax.annotate(f'{np.around(time_snapshot, DECIMAL_COUNT)} Gyr', xy=(-width/YEAR_X_FRAC, width/YEAR_Y_FRAC), color='black', zorder=4, fontsize=DECIMAL_FONT_SIZE)

    if SHOW_CBAR:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3.5%", pad=0.1)
        cb = plt.colorbar(im, cax=cax)
        cb.set_label(r'$\log \, \rho \;$ (g cm$^{-3}) $', labelpad=8)

    plt.savefig(f'{IMAGE_PATH}density_{count:03d}.png', dpi=300)
    plt.close()

    # DM DENSITY ---------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1)

    plt.sca(ax)  # Define o eixo atual para o sph.image
    im_array = sph.image(
        df.dm,
        qty="rho",
        units="g cm^-3",
        width=width,
        cmap="twilight",
        av_z='rho',
        vmin=2e-28, 
        vmax=3e-25,
        show_cbar=SHOW_CBAR,
        noplot=True
    )

    # Constrói o plot manualmente
    extent = (-width/2, width/2, -width/2, width/2)
    norm = clrs.LogNorm(vmin=1e-28, vmax=1e-24)
    im = ax.imshow(im_array, cmap="twilight", extent=extent, norm=norm)
    ax.invert_xaxis()

    fig.subplots_adjust(left=LEFT, bottom=BOTTOM, top=TOP, right=RIGHT, hspace=HSPACE, wspace=WSPACE)
    ax.set_xlabel(r'$x$ (kpc)')
    ax.set_ylabel(r'$y$ (kpc)')
    ax.set_aspect('equal')
    ax.annotate(f'{np.around(time_snapshot, DECIMAL_COUNT)} Gyr', xy=(-width/YEAR_X_FRAC, width/YEAR_Y_FRAC), color='black', zorder=4, fontsize=DECIMAL_FONT_SIZE)

    if SHOW_CBAR:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3.5%", pad=0.1)
        cb = plt.colorbar(im, cax=cax)
        cb.set_label(r'$\log \, \rho \;$ (g cm$^{-3}) $', labelpad=8)

    plt.savefig(f'{IMAGE_PATH}dm_{count:03d}.png', dpi=300)
    plt.close()

    # TEMPERATURE ----------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1)

    plt.sca(ax)  # Define o eixo atual para o sph.image
    im_array = sph.image(
        df.gas,
        qty="kT",
        width=width,
        cmap="inferno",
        show_cbar=SHOW_CBAR,
        noplot=True
    )

    # Constrói o plot manualmente
    extent = (-width/2, width/2, -width/2, width/2)
    norm = clrs.Normalize(vmin=0, vmax=12)
    # norm = clrs.Normalize(vmin=0, vmax=15)
    im = ax.imshow(im_array, cmap="inferno", extent=extent, norm=norm)
    ax.invert_xaxis()

    fig.subplots_adjust(left=LEFT, bottom=BOTTOM, top=TOP, right=RIGHT, hspace=HSPACE, wspace=WSPACE)
    ax.set_xlabel(r'$x$ (kpc)')
    ax.set_ylabel(r'$y$ (kpc)')
    ax.set_aspect('equal')
    ax.annotate(f'{np.around(time_snapshot, DECIMAL_COUNT)} Gyr', xy=(-width/YEAR_X_FRAC, width/YEAR_Y_FRAC), color='white', zorder=4, fontsize=DECIMAL_FONT_SIZE)

    if SHOW_CBAR:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3.5%", pad=0.1)
        cb = plt.colorbar(im, cax=cax)
        cb.set_label(r'$kT$ (keV)', labelpad=8)

    plt.savefig(f'{IMAGE_PATH}temperature_{count:03d}.png', dpi=300)
    plt.close()

    count += 1

# CRIAR ANIMAÇÕES ----------------------------------------------------------------------

# ffmpeg -framerate 12 -i dm_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y dm.mp4
# ffmpeg -framerate 12 -i density_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y density.mp4
# ffmpeg -framerate 12 -i temperature_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y temperature.mp4

# ffmpeg -i dm.mp4 -i density.mp4 -i temperature.mp4 -filter_complex "hstack=inputs=3" -c:v libx264 -pix_fmt yuv420p -y combined.mp4

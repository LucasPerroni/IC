import math
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PATH_TO_HDF5 = "./hdf5/"
PATH_TO_PLOT = "./plots/"
PATH_TO_JOIN_FRAMES = "./frames/join/"
lines = open(PATH_TO_HDF5 + "snapshot.txt", "r").readlines()
interactions = 20

mi = 0.6 # Average molecular weight
Mh = 1.67262192 * 10**(-27) # Proton mass in kg

# Cluster halo variables
Mh_cluster = 50000 # halo mass
Ah = 200 # scale factor in the Dehnen density profile
Gamma_h = 1 # Dehnen density profile free parameter. If gamma = 1, then the profile
            # is equal to a Hernquist profile, featuring a central cusp. If gamma = 0,
            # then the profile features a central core.

# Cluster gas variables
Mg = 15000 # gas mass
Ag = 200 # scale factor in the Dehnen density profile
Gamma_g = 0

# Analytic function data
min_radius = 0.1 # first point of the array
step = 0.1

# cores
pontos = '#4772FF'
linha = '#BACBFF'

# Carrega uma imagem para obter suas dimensões
img = mpimg.imread(f"{PATH_TO_JOIN_FRAMES}out_00001.jpg")
img_height, img_width, _ = img.shape  # Altura e largura em pixels
dpi = 100
fig_width = img_width / dpi
fig_height = img_height / dpi

for i in range(len(lines)):
    print(f"Calculating snapshot_{i:03d}")
    density = []
    density_c= []
    radius = []
    radius_c = []

    temperature = []
    density2 = []
    density2_c = []
    radius2 = []
    radius2_c = []
    radius2_log = []

    snapshot = f'{PATH_TO_HDF5}snapshot_{i:03d}.hdf5'
    s = h5py.File(snapshot, 'r')

    t = s['Header'].attrs[u'Time']

    mass = s['PartType1']['Masses'][:] * 1e10 # Halo
    x = s['PartType1']['Coordinates'][:,0] # Halo
    y = s['PartType1']['Coordinates'][:,1] # Halo
    z = s['PartType1']['Coordinates'][:,2] # Halo

    mass2 = s['PartType0']['Masses'][:] * 1e10 # Gas
    x2 = s['PartType0']['Coordinates'][:,0] # Gas
    y2 = s['PartType0']['Coordinates'][:,1] # Gas
    z2 = s['PartType0']['Coordinates'][:,2] # Gas
    u_tot = s['PartType0']['InternalEnergy'][:] * 10**6 # (km/s)² -> (m/s)²

    # Halo
    R = np.sqrt(x**2 + y**2 + z**2)
    Rmax = np.max(R)

    # Gas
    R2 = np.sqrt(x2**2 + y2**2 + z2**2)
    Rmax2 = np.max(R2)

    # Define the number of bins for the logarithmic radial bins
    n_bins = interactions + 1

    # logarithmically spaced radial bins
    log_bins = np.logspace(np.log10(1), np.log10(Rmax), n_bins)
    log_g_bins = np.logspace(np.log10(1), np.log10(Rmax2), n_bins)

    # density and radius for gas and halo in logarithmic scale
    for j in range(n_bins - 1):
        rmin = log_bins[j]
        rmax = log_bins[j + 1]

        rmin2 = log_g_bins[j]
        rmax2 = log_g_bins[j + 1]

        # Halo
        cond = np.argwhere(np.logical_and(R > rmin, R < rmax))
        m = np.sum(mass[cond])
        vol = (4/3 * math.pi * rmax**3) - (4/3 * math.pi * rmin**3)
        d = m / vol
        r = (rmax + rmin) / 2
        density.append(d)
        radius.append(r)

        # Gas
        cond2 = np.argwhere(np.logical_and(R2 > rmin2, R2 < rmax2))
        m2 = np.sum(mass2[cond2])
        vol2 = (4/3 * math.pi * rmax2**3) - (4/3 * math.pi * rmin2**3)
        d2 = m2 / vol2
        r2 = (rmax2 + rmin2) / 2
        density2.append(d2)
        radius2_log.append(r2)

    # temperature of gas
    passo2 = Rmax2/interactions
    rmin2 = 0
    rmax2 = passo2
    for k in range(interactions):
        cond2 = np.argwhere(np.logical_and(R2 > rmin2, R2 < rmax2))
        u = np.sum(u_tot[cond2]) / len(u_tot[cond2]) # Average energy density
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV
        r2 = (rmax2 + rmin2) / 2
        temperature.append(kT)
        radius2.append(r2)
        rmin2 += passo2
        rmax2 += passo2
    
    # Analytic functions for density and radius
    for r in np.arange(min_radius, Rmax, step):
        p = ((3 - Gamma_h) * Mh_cluster * 1e10 / (4 * math.pi)) * (Ah / (r**Gamma_h * (r + Ah)**(4 - Gamma_h)))
        density_c.append(p)
        radius_c.append(r)
    for r2 in np.arange(min_radius, Rmax2, step):
        p2 = ((3 - Gamma_g) * Mg * 1e10 / (4 * math.pi)) * (Ag / (r2**Gamma_g * (r2 + Ag)**(4 - Gamma_g)))
        density2_c.append(p2)
        radius2_c.append(r2)

    # plot
    if Gamma_h == 0:
        funcao = "Dehnen density profile"
    elif Gamma_h == 1:
        funcao = "Hernquist density profile"

    if Gamma_g == 0:
        funcao2 = "Dehnen density profile"
    elif Gamma_g == 1:
        funcao2 = "Hernquist density profile"

    fig = plt.figure(figsize=(fig_width, fig_height * 5/3))
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    ax1.plot(radius, density, '.', ms = 12, mec = pontos, mfc = pontos, label='Density from Snapshot')
    ax1.plot(radius_c, density_c, '-', ms = 1.5, color = linha, label=funcao)
    ax1.set_title(f'(Halo) Density x Radius - {t} Gyr')
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel('$\\rho$ ($M_{\odot} kpc^{-3}$)')
    ax1.set_xlabel('$r$ ($kpc$)')
    ax1.set_xlim(1e0, 4e3)
    ax1.set_ylim(1e1, 1e9)
    ax1.set_aspect('auto')
    ax1.legend()

    ax2.plot(radius2_log, density2, '.', ms = 12, mec = pontos, mfc = pontos, label='Density from Snapshot')
    ax2.plot(radius2_c, density2_c, '-', ms = 1.5, color = linha, label=funcao2)
    ax2.set_title(f'(Gas) Density x Radius - {t} Gyr')
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_ylabel('$\\rho$ ($M_{\odot} kpc^{-3}$)')
    ax2.set_xlabel('$r$ ($kpc$)')
    ax2.set_xlim(1e0, 4e3)
    ax2.set_ylim(1e2, 1e7)
    ax2.set_aspect('auto')
    ax2.legend()

    ax3.plot(radius2, temperature, '.', ms = 12, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    ax3.set_title(f'Temperature x Radius - {t} Gyr')
    ax3.set_ylabel('$kT$ ($keV$)')
    ax3.set_xlabel('$r$ ($kpc$)')
    ax3.set_xlim(0, 3000)
    ax3.set_ylim(0, 12)
    ax3.set_aspect('auto')
    ax3.legend()

    plt.savefig(f'{PATH_TO_PLOT}d-r_{i:03d}.jpg')
    plt.close(fig)

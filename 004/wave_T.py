import math
import h5py
import numpy as np
import matplotlib.pyplot as plt

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

SNAPSHOT_PATH = "/mnt/d/IC/snapshots/"

mi = 0.6 # Average molecular weight
Mh = 1.67262192 * 10**(-27) # Proton mass in kg

# cores
pontos = '#4772FF'
linha = '#BACBFF'

lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    # if i < 7 or i > 12:
    #     continue

    print(f"{bcolors.OKGREEN}Plotting snapshot {i:03d} of {len(lines) - 1:03d}...{bcolors.ENDC}")
    snapshot = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    s = h5py.File(snapshot, 'r')

    time = s['Header'].attrs[u'Time'] # Gyr
    mass = s['PartType0']['Masses'][:] # 1e10 Msol
    u_tot = s['PartType0']['InternalEnergy'][:] * 10**6 # (km/s)² -> (m/s)²

    x = s['PartType0']['Coordinates'][:,0] # kpc
    y = s['PartType0']['Coordinates'][:,1] # kpc
    z = s['PartType0']['Coordinates'][:,2] # kpc
    R = np.sqrt(x**2 + y**2 + z**2) # radius

    x_plot = []
    kT_plot = []

    interacoes = 200
    length = max(x) - min(x)
    for j in range(interacoes):
        x1 = min(x) + (j * length / interacoes)
        x2 = min(x) + ((j + 1) * length / interacoes)
        cond = np.argwhere((x > x1) & (x < x2) & (y > -200) & (y < 200) & (z > -200) & (z < 200))

        u = np.sum(u_tot[cond]) / len(u_tot[cond])
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

        x_plot.append((x1 + x2) / 2)
        kT_plot.append(kT)

    # interacoes = 200
    # length = 1000 - 300
    # for j in range(interacoes):
    #     x1 = 300 + (j * length / interacoes)
    #     x2 = 300 + ((j + 1) * length / interacoes)
    #     cond = np.argwhere((x > x1) & (x < x2) & (y > -200) & (y < 200) & (z > -200) & (z < 200))

    #     u = np.sum(u_tot[cond]) / len(u_tot[cond])
    #     kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

    #     x_plot.append((x1 + x2) / 2)
    #     kT_plot.append(kT)


    # Plot
    fig, ax1 = plt.subplots(figsize=(8, 6))

    # ax1 - temperature x radius
    ax1.plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    ax1.set_title(f"Temperature x Radius - {time} Gyr")
    ax1.set_ylabel('$kT$ ($keV$)')
    ax1.set_xlabel('$x$ ($kpc$)')
    ax1.set_xlim(-3000, 3000)
    # ax1.set_xlim(300, 900)
    ax1.set_ylim(0, 40)
    ax1.set_aspect('auto')
    ax1.legend()

    # save figure
    plt.savefig(f"plot/t-r_{i:03d}.png")
    # plt.savefig(f"plot/wave/t-r_{i:03d}.png")

# ffmpeg -framerate 6 -i plot/t-r_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

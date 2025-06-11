import h5py
import numpy as np
from sympy import *
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

# parâmetros para plot de 4 imagens
plot_index = 0
fig, axs = plt.subplots(2, 2, figsize=(16, 12))

lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    if i < 8 or i > 11:
    # if i != 8:
        continue

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

    # interacoes = 200
    # length = max(x) - min(x)
    # for j in range(interacoes):
    #     x1 = min(x) + (j * length / interacoes)
    #     x2 = min(x) + ((j + 1) * length / interacoes)
    #     cond = np.argwhere((x > x1) & (x < x2) & (y > -200) & (y < 200) & (z > -200) & (z < 200))

    #     u = np.sum(u_tot[cond]) / len(u_tot[cond])
    #     kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

    #     x_plot.append((x1 + x2) / 2)
    #     kT_plot.append(kT)

    interacoes = 200
    length = 1000 - 300
    for j in range(interacoes):
        x1 = 300 + (j * length / interacoes)
        x2 = 300 + ((j + 1) * length / interacoes)
        cond = np.argwhere((x > x1) & (x < x2) & (y > -100) & (y < 100) & (z > -100) & (z < 100))

        u = np.sum(u_tot[cond]) / len(u_tot[cond])
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

        x_plot.append((x1 + x2) / 2)
        kT_plot.append(kT)

    # Medindo a descontinuidade de temperatura
    sep = 8 # numero de pontos entre T1/T2 e a descontinuidade 
    grad_kT = np.gradient(kT_plot, x_plot)
    idx = np.where(grad_kT == np.min(grad_kT))[0][0]
    idx_menor = idx + sep
    idx_maior = idx - sep
    T1 = kT_plot[idx_menor]
    T2 = kT_plot[idx_maior]

    # Calculando número de Mach
    M = symbols('M')
    eq = Eq((5*M**4 + 14*M**2 - 3) / (16*M**2), T2/T1)
    mach = solve(eq)
    mach = max([m for m in mach if im(m) == 0])

    # Plot
    # fig, ax = plt.subplots(figsize=(8, 6))
    # ax.plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    # ax.set_title(f"Temperature x Radius - {time} Gyr")
    # ax.set_ylabel('$kT$ ($keV$)')
    # ax.set_xlabel('$x$ ($kpc$)')
    # # ax.set_xlim(-3000, 3000)
    # ax.set_xlim(300, 900)
    # ax.set_ylim(0, 40)
    # ax.set_aspect('auto')
    # ax.axvline(x=x_plot[idx], color='r', linestyle='--', alpha=0.5, label=f'Descontinuidade')
    # ax.legend()

    i_plot = plot_index // 2
    j_plot = plot_index % 2
    axs[i_plot][j_plot].plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    axs[i_plot][j_plot].set_title(f"Temperature x Radius - {time:.3f} Gyr - Mach: {mach:.4f}")
    axs[i_plot][j_plot].set_ylabel('$kT$ ($keV$)')
    axs[i_plot][j_plot].set_xlabel('$x$ ($kpc$)')
    axs[i_plot][j_plot].set_xlim(300, 900)
    axs[i_plot][j_plot].set_ylim(0, 40)
    axs[i_plot][j_plot].set_aspect('auto')
    axs[i_plot][j_plot].axvline(x=x_plot[idx], color='r', linestyle='--', alpha=0.5, label=f'Descontinuidade')
    axs[i_plot][j_plot].axhline(y=kT_plot[idx_maior], color='g', linestyle='--', alpha=0.5, label=f'T2')
    axs[i_plot][j_plot].axhline(y=kT_plot[idx_menor], color='purple', linestyle='--', alpha=0.5, label=f'T1')
    # axs[i_plot][j_plot].axvline(x=x_plot[idx_maior], color='g', linestyle='--', alpha=0.5, label=f'T2')
    # axs[i_plot][j_plot].axvline(x=x_plot[idx_menor], color='purple', linestyle='--', alpha=0.5, label=f'T1')
    axs[i_plot][j_plot].legend()
    plot_index += 1
    
    # save figure
    # plt.savefig(f"plot/t-r_{i:03d}.png")

plt.savefig(f"plot/t-r_all.png")
# ffmpeg -framerate 6 -i plot/t-r_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

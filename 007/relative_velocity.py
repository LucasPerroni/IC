import sys
import h5py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.stats import binned_statistic_2d

def densidade_pico(x_coords, y_coords, bins=1000):
    H, xedges, yedges, binnumber = binned_statistic_2d(
        x_coords, y_coords, None, statistic='count', bins=bins
    )

    idx_max = np.unravel_index(np.argmax(H), H.shape)
    x_center = 0.5 * (xedges[idx_max[0]] + xedges[idx_max[0]+1])
    y_center = 0.5 * (yedges[idx_max[1]] + yedges[idx_max[1]+1])
    return x_center, y_center

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

SNAPSHOT_CODE = "0005_0006_100/"
SNAPSHOT_PATH = "/mnt/d/UFPR/IC/snapshots/" + SNAPSHOT_CODE
IMAGE_PATH = "plots/" + SNAPSHOT_CODE
MACH_ANALISYS = True
FONT_SIZE = 14

# cores
pontos = '#4772FF'
linha = '#BACBFF'

cluster_1 = {"x": [], "y": []}
cluster_2 = {"x": [], "y": []}
distance_between = {"dist": []}
analised_timestamp = {"time_init": None, "time_end": None, "dist_init": None, "dist_end": None}
sim_time = []

dist_prev = None
shock_detected = False
if MACH_ANALISYS:
    s_init, s_end = 90, 100 # 0005_0006_100
    # s_init, s_end = 60, 70 # 0005_0006_1000
    # s_init, s_end = 40, 50 # 0005_0006_2000

lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    # if i % 3 != 0:
    #     continue

    print(f"{bcolors.OKGREEN}Plotting snapshot {i:03d} of {len(lines) - 1:03d}...{bcolors.ENDC}")
    snapshot = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    s = h5py.File(snapshot, 'r')

    time = s['Header'].attrs[u'Time'] # Gyr
    mass = s['PartType1']['Masses'][:] # 1e10 Msol

    x = s['PartType1']['Coordinates'][:,0] # kpc
    x = -x
    y = s['PartType1']['Coordinates'][:,1] # kpc
    z = s['PartType1']['Coordinates'][:,2] # kpc

    # Calculando as coordenadas x e y de cada pico de densidade
    x_median = np.median(x)
    mask_halo1 = x < x_median
    mask_halo2 = x >= x_median

    x1, y1 = x[mask_halo1], y[mask_halo1]
    x2, y2 = x[mask_halo2], y[mask_halo2]

    cx1, cy1 = densidade_pico(x1, y1)
    cx2, cy2 = densidade_pico(x2, y2)
    dist = np.sqrt((cx2 - cx1)**2 + (cy2 - cy1)**2)  # kpc

    # Inverte as coordenadas após o choque
    if (dist_prev is not None) & (not shock_detected):
        if dist > dist_prev + 20:
            # print("foi")
            shock_detected = True
            cx1, cx2, cy1, cy2 = cx2, cx1, cy2, cy1
    elif shock_detected:
        cx1, cx2, cy1, cy2 = cx2, cx1, cy2, cy1
    dist_prev = dist

    cluster_1['x'].append(cx1)
    cluster_2['x'].append(cx2)
    distance_between['dist'].append(dist)
    sim_time.append(time)

    # Pega o tempo e a distância de separação no momento de início e fim da análise do número de Mach
    if MACH_ANALISYS & (i == s_init):
        analised_timestamp['time_init'] = time
        analised_timestamp['dist_init'] = dist
    elif MACH_ANALISYS & (i == s_end):
        analised_timestamp['time_end'] = time
        analised_timestamp['dist_end'] = dist

# ---- Plots Cluster 1
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(sim_time, cluster_1['x'], '.', ms = 4, mec = pontos, mfc = pontos, label='Centro do cluster')
# ax.set_title(f"Centro do Cluster 1 x Tempo")
ax.set_ylabel('$x$ (kpc)', fontsize=FONT_SIZE)
ax.set_xlabel('$t$ (Gyr)', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
# ax.legend()
plt.savefig(f"{IMAGE_PATH}pos-time-1.png")

# ---- Plots Cluster 2
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(sim_time, cluster_2['x'], '.', ms = 4, mec = pontos, mfc = pontos, label='Centro do cluster')
# ax.set_title(f"Centro do Cluster 2 x Tempo")
ax.set_ylabel('$x$ (kpc)', fontsize=FONT_SIZE)
ax.set_xlabel('$t$ (Gyr)', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
# ax.legend()
plt.savefig(f"{IMAGE_PATH}pos-time-2.png")

# ---- Plot Distância entre Clusters
dist_smooth = savgol_filter(distance_between['dist'], window_length=9, polyorder=3)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(sim_time, dist_smooth, '.', ms = 4, mec = pontos, mfc = pontos, label='Distância absoluta')
# ax.set_title(f"Distância entre os centros")
ax.set_ylabel('$x$ (kpc)', fontsize=FONT_SIZE)
ax.set_xlabel('$t$ (Gyr)', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
if MACH_ANALISYS:
    ax.axvline(analised_timestamp['time_init'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvline(analised_timestamp['time_end'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvspan(analised_timestamp['time_init'], analised_timestamp['time_end'], facecolor='#aba9a9', alpha=0.2, label="Intervalo de análise do Número de Mach")
# ax.legend(loc="upper left", fontsize=FONT_SIZE)
plt.tight_layout()
plt.savefig(f"{IMAGE_PATH}dist.png")

# ---- Plot Velocidade x Tempo
velocity = np.gradient(dist_smooth, sim_time)*0.9778 # kpc/gyr -> km/s
cond = (sim_time >= analised_timestamp['time_init']) & (sim_time <= analised_timestamp['time_end'])
analised_v = np.mean(velocity[cond])

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(sim_time, velocity, '.', ms = 4, mec = pontos, mfc = pontos, label='Velocidade relativa')
# ax.set_title(f"Velocidade Relativa x Tempo")
ax.set_ylabel('$vx$ (km/s)', fontsize=FONT_SIZE)
ax.set_xlabel('$t$ (Gyr)', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
ax.axhline(0, linestyle='--', color="#999797", alpha=0.2)
ax.axhline(analised_v, ls='-', color="#363636", alpha=0.2, label="Velocidade relativa média na análise de Mach")
if MACH_ANALISYS:
    ax.axvline(analised_timestamp['time_init'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvline(analised_timestamp['time_end'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvspan(analised_timestamp['time_init'], analised_timestamp['time_end'], facecolor='#aba9a9', alpha=0.2, label="Intervalo de análise do Número de Mach")
# ax.legend(loc="lower right", fontsize=FONT_SIZE)
plt.tight_layout()
plt.savefig(f"{IMAGE_PATH}vel-time.png")

# ---- Plot Velocidade x Distância
dist_signed = np.sign(np.array(cluster_1['x']) - np.array(cluster_2['x'])) * dist_smooth
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(dist_signed, velocity, '.', ms = 4, mec = pontos, mfc = pontos, label='Velocidade relativa')
# ax.set_title(f"Velocidade Relativa x Distância entre os centros")
ax.set_ylabel('$vx$ (km/s)', fontsize=FONT_SIZE)
ax.set_xlabel('$x$ (kpc)', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
ax.axhline(0, linestyle='--', color="#999797", alpha=0.2)
ax.axhline(analised_v, ls='-', color="#363636", alpha=0.2, label="Velocidade relativa média na análise de Mach")
if MACH_ANALISYS:
    ax.axvline(analised_timestamp['dist_init'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvline(analised_timestamp['dist_end'], linestyle='--', color="#aba9a9", alpha=0.2)
    ax.axvspan(analised_timestamp['dist_init'], analised_timestamp['dist_end'], facecolor='#aba9a9', alpha=0.2, label="Intervalo de análise do Número de Mach")
# ax.legend(loc="lower right", fontsize=FONT_SIZE)
plt.tight_layout()
plt.savefig(f"{IMAGE_PATH}vel-pos.png")

print(f"Velocidade média: {analised_v:.4f}")

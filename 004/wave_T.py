import sys
import h5py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

def velocidade_som_keV(kT_keV, mi=0.6, Mh=1.67262192e-27, gamma=5/3):
    """
    Calcula a velocidade do som a partir de kT em keV.
    """
    kT_joule = kT_keV * 1.60218e-16  # keV -> J
    cs = np.sqrt(gamma * kT_joule / (mi * Mh))  # m/s
    return cs / 1000  # em km/s

def fill_nan_nearest(arr):
    arr = np.array(arr, dtype=float)
    nans = np.isnan(arr)
    if not np.any(nans):
        return arr

    idx = np.arange(len(arr))
    valid_idx = idx[~nans]
    valid_values = arr[~nans]

    nearest_idx = np.abs(valid_idx[:, None] - idx).argmin(axis=0)
    arr[nans] = valid_values[nearest_idx[nans]]
    return arr


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
IMAGE_PATH = "plot/" + SNAPSHOT_CODE
PLOT_INFOS = True
FULL_PLOT = not PLOT_INFOS
FONT_SIZE = 14

mi = 0.6 # Average molecular weight
Mh = 1.67262192 * 10**(-27) # Proton mass in kg
gamma_var = 5/3 # coeficiente de dilatação adiabática do gás

# cores
pontos = '#4772FF'
linha = '#BACBFF'

dxdt = {"x": [], "t": []} # objeto para plotar a posição da descontinuidade pelo tempo
velocities = {"cs": [], "u": []} # objeto para calcular o número de mach pela velocidade
machs = {"time": [], "mach": []} # objeto para guardar os machs calculados por T2/T1
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
s_init, s_end = 90, 100 # 0005_0006_100
# s_init, s_end = 60, 70 # 0005_0006_1000
# s_init, s_end = 40, 50 # 0005_0006_2000
for i in range(len(lines)):
    if (i < s_init) | (i > s_end):
        continue

    print(f"{bcolors.OKGREEN}Plotting snapshot {i:03d} of {len(lines) - 1:03d}...{bcolors.ENDC}")
    snapshot = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    s = h5py.File(snapshot, 'r')

    time = s['Header'].attrs[u'Time'] # Gyr
    mass = s['PartType0']['Masses'][:] # 1e10 Msol
    u_tot = s['PartType0']['InternalEnergy'][:] * 10**6 # (km/s)² -> (m/s)²

    vx = s['PartType0']['Velocities'][:,0]
    vx = -vx
    vy = s['PartType0']['Velocities'][:,1]
    vz = s['PartType0']['Velocities'][:,2]
    V = np.sqrt(vx**2 + vy**2 + vz**2)

    x = s['PartType0']['Coordinates'][:,0] # kpc
    x = -x
    y = s['PartType0']['Coordinates'][:,1] # kpc
    z = s['PartType0']['Coordinates'][:,2] # kpc
    R = np.sqrt(x**2 + y**2 + z**2) # radius

    x_plot = []
    kT_plot = []
    limit_yz = (y > -100) & (y < 100) & (z > -100) & (z < 100)

    # Calcular kT para o plot de temperatura
    interacoes = 145 # 0005_0006_100
    # interacoes = 178 # 0005_0006_1000
    # interacoes = 110 # 0005_0006_2000
    if FULL_PLOT:
        lim_sup = 1000
        lim_inf = -1000
        length = lim_sup - lim_inf
        for j in range(interacoes):
            x1 = lim_inf + (j * length / interacoes)
            x2 = lim_inf + ((j + 1) * length / interacoes)
            cond = (x > x1) & (x < x2) & limit_yz

            u = np.mean(u_tot[cond])
            kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

            x_plot.append((x1 + x2) / 2)
            kT_plot.append(kT)
    else:
        length = 1000 - 200
        for j in range(interacoes):
            x1 = 200 + (j * length / interacoes)
            x2 = 200 + ((j + 1) * length / interacoes)
            cond = (x > x1) & (x < x2) & limit_yz

            u = np.mean(u_tot[cond])
            kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV

            x_plot.append((x1 + x2) / 2)
            kT_plot.append(kT)
    kT_plot = fill_nan_nearest(kT_plot)

    if PLOT_INFOS:
        # Medindo a descontinuidade de temperatura
        sep = 18 # numero de pontos entre T1/T2 e a descontinuidade

        x_plot = np.array(x_plot)
        mask_range = (x_plot > 300) & (x_plot < 660) # 0005_0006_100
        # mask_range = (x_plot > 300) & (x_plot < 840) # Resto
        grad_kT = np.full_like(kT_plot, 0, dtype=float)
        grad_kT[mask_range] = np.gradient(kT_plot[mask_range], x_plot[mask_range])

        idx = np.where(grad_kT == np.min(grad_kT))[0][0]
        idx_menor = idx + sep
        idx_maior = np.where(kT_plot == np.max(kT_plot))[0][0]

        T1 = kT_plot[idx_menor]
        T2 = kT_plot[idx_maior]

        # Calculando as velocidades no gás não chocado
        plot_limit = limit_yz & (x > 200)
        descontinuity_limit = limit_yz & (x > x_plot[idx])
        # x_nao_chocado = sorted(x[descontinuity_limit])[int(len(x[descontinuity_limit])/1.05)]
        x_nao_chocado = x_plot[idx] + 250
        volume_limit = limit_yz & (x > x_nao_chocado -100) & (x < x_nao_chocado + 100)

        u = np.mean(u_tot[volume_limit]) # Energia interna no gás não chocado
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV
        cs = velocidade_som_keV(kT) # velocidade do som no gás não chocado
        u = np.mean(vx[volume_limit]) # velocidade do gás não chocado

        velocities['cs'].append(cs)
        velocities['u'].append(u)

        # fig, ax = plt.subplots(figsize=(8, 6))
        # ax.plot(x[plot_limit], vx[plot_limit], ".", markersize=1.2)
        # # ax.set_title(f"Velocity x Radius - {time:.3f} Gyr")
        # ax.set_ylabel("velocity (km/s)")
        # ax.set_xlabel("x (kpc)")
        # ax.axvline(x=x_nao_chocado, color="y", linestyle="--", alpha=0.5, label="Não Chocado")
        # ax.axvline(x=x_plot[idx], color="black", linestyle="--", alpha=0.5, label="Descontinuidade")
        # ax.axhline(y=0, color='r', linestyle='--', alpha=0.5, label="0 km/s")
        # ax.set_xlim(200, 1200)
        # ax.set_ylim(-600, 1800)
        # ax.legend()
        # plt.tight_layout()
        # plt.savefig(f"{IMAGE_PATH}x_vx_{i:03d}.png")

        # Calculando número de Mach
        if T1 == 0.0:
            T1 = 0.1
        M = symbols('M')
        eq = Eq((5*M**4 + 14*M**2 - 3) / (16*M**2), T2/T1)
        mach = solve(eq)
        mach = max([m for m in mach if im(m) == 0])

        machs['time'].append(time)
        machs['mach'].append(mach)

        dxdt["x"].append(x_plot[idx])
        dxdt['t'].append(time)

    # Plot individual
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperatura das partículas')
    # ax.set_title(f"Temperature x Radius - {time:.3f} Gyr")
    ax.set_ylabel('$kT$ (keV)', fontsize=FONT_SIZE)
    ax.set_xlabel('$x$ (kpc)', fontsize=FONT_SIZE)
    ax.tick_params(axis='both', labelsize=FONT_SIZE)
    if FULL_PLOT:
        ax.set_xlim(lim_inf, lim_sup)
        ax.axhline(y=max(kT_plot), color='g', linestyle='--', alpha=0.5, 
                   label=f'Maior temperatura: {max(kT_plot):.2f}')
    else:
        ax.set_xlim(200, 1000)
        ax.set_ylim(0, 10)
    if PLOT_INFOS:
        # ax.axhline(y=kT_plot[idx_maior], color='g', linestyle='--', alpha=0.2, label=f'T2')
        # ax.axhline(y=kT_plot[idx_menor], color='purple', linestyle='--', alpha=0.2, label=f'T1')
        ax.axvline(x=x_plot[idx], ls='-', color="#363636", alpha=0.2, label=f'Descontinuidade')
    ax.set_aspect('auto')
    # ax.legend(loc="upper right", fontsize=FONT_SIZE)
    plt.text(
        0.98, 0.98,           # posição (x, y) em coordenadas relativas ao eixo
        f"{time:.3f} Gyr",       # texto
        transform=plt.gca().transAxes,  # garante que a posição é relativa ao gráfico
        ha="right", va="top", # alinhamento
        color="#363636", fontsize=FONT_SIZE, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none") # fundo branco sem borda
    )
    plt.tight_layout()
    if FULL_PLOT:
        plt.savefig(f"{IMAGE_PATH}full_plot/t-r_{i:03d}.png")
    else:
        plt.savefig(f"{IMAGE_PATH}t-r_{i:03d}.png")

if PLOT_INFOS:
    # Ajuste linear a curva de posição pelo tempo
    a, b = np.polyfit(dxdt['t'], dxdt['x'], 1)
    x_ajuste = np.linspace(min(dxdt['t']), max(dxdt['t']), 100)
    y_ajuste = a*x_ajuste + b

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(dxdt['t'], dxdt["x"], ".", label="Posição da descontinuidade")
    ax.plot(x_ajuste, y_ajuste, color="lightblue", label="Ajuste linear")
    ax.set_xlabel("$t$ (Gyr)", fontsize=FONT_SIZE)
    ax.set_ylabel("$x$ (kpc)", fontsize=FONT_SIZE)
    ax.tick_params(axis='both', labelsize=FONT_SIZE)
    # ax.set_title(f"Posição da Descontinuidade x Tempo - Velocidade: {a*0.9778:.2f} km/s")
    # ax.legend(loc="lower right", fontsize=FONT_SIZE)
    plt.tight_layout()
    plt.savefig(f"{IMAGE_PATH}pos-time.png")

    # Print dos outputs
    print(f"Velocidade da onda: {a*0.9778:.2f} km/s")
    print(f"Velocidade do som não chocado: {np.mean(velocities['cs']):.2f} km/s")
    print(f"Velocidade das particulas não chocadas: {np.mean(velocities['u']):.2f} km/s")

    mach_v = (a*0.9778 - np.mean(velocities['u'])) / np.mean(velocities['cs'])
    print(f"Mach pela velocidade: {mach_v:.2f} \n")

    for i in range(len(machs["time"])):
        print(f"Mach por T2/T1 em {machs['time'][i]:.3f} Gyr: {machs['mach'][i]:.2f}")
    print(f"Média dos machs por T2/T1: {np.mean(machs['mach']):.2f}")

# ffmpeg -framerate 6 -i plot/t-r_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

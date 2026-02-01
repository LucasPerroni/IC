import os
import h5py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from functions import *

SNAPSHOT_CODE = "0005_0006_0/"
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
erro_savgol_mach = {"time": [], "erro": []} # objeto para guardar os erros dos machs pelo filtro Savitzky-Golay
last_snapshot_descontinuity = 0 # pega a posição da descontinuidade no ultimo snapshot
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
os.makedirs(f"{IMAGE_PATH}", exist_ok=True)

# CONFIGURAÇÕES POR SNAPSHOT
SNAPSHOT_CONFIG = {
    "0005_0005_0/": {
        "s_init": 80,
        "s_end": 90,
        "mask_lower": -100,
        "mask_upper": 400,
        "interacoes": 200,
        "x_init": -200,
        "x_end": 600,
        "max_kT": 17,
        "max_dx": 100,
    },
    "0005_0005_1000/": {
        "s_init": 54,
        "s_end": 64,
        "mask_lower": 400,
        "mask_upper": 1050,
        "interacoes": 200,
        "x_init": 300,
        "x_end": 1200,
        "max_kT": 20,
        "max_dx": 100,
    },
    "0005_0005_3000/": {
        "s_init": 30,
        "s_end": 40,
        "mask_lower": 600,
        "mask_upper": 1700,
        "interacoes": 200,
        "x_init": 600,
        "x_end": 1900,
        "max_kT": 30,
        "max_dx": 200,
    },
    "0005_0006_0/": {
        "s_init": 95,
        "s_end": 105,
        "mask_lower": 270,
        "mask_upper": 660,
        "interacoes": 200,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 12,
        "max_dx": 100,
    },
    "0005_0006_100/": {
        "s_init": 90,
        "s_end": 100,
        "mask_lower": 270,
        "mask_upper": 660,
        "interacoes": 200,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 10,
        "max_dx": 100,
    },
    "0005_0006_1000/": {
        "s_init": 60,
        "s_end": 70,
        "mask_lower": 270,
        "mask_upper": 840,
        "interacoes": 200,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 10,
        "max_dx": 100,
    },
    "0005_0006_2000/": {
        "s_init": 40,
        "s_end": 50,
        "mask_lower": 270,
        "mask_upper": 840,
        "interacoes": 200,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 10,
        "max_dx": 100,
    },
    "0005_0006_3000/": {
        "s_init": 30,
        "s_end": 40,
        "mask_lower": 270,
        "mask_upper": 840,
        "interacoes": 150,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 12,
        "max_dx": 100,
    },
    "0005_0007_0/": {
        "s_init": 82,
        "s_end": 92,
        "mask_lower": -50,
        "mask_upper": 400,
        "interacoes": 200,
        "x_init": -100,
        "x_end": 600,
        "max_kT": 15,
        "max_dx": 100,
    },
    "0005_0007_1000/": {
        "s_init": 53,
        "s_end": 63,
        "mask_lower": 270,
        "mask_upper": 900,
        "interacoes": 200,
        "x_init": 200,
        "x_end": 1000,
        "max_kT": 17,
        "max_dx": 100,
    },
    "0005_0007_3000/": {
        "s_init": 30,
        "s_end": 40,
        "mask_lower": 500,
        "mask_upper": 1500,
        "interacoes": 200,
        "x_init": 450,
        "x_end": 1600,
        "max_kT": 25,
        "max_dx": 200,
    },
}

try:
    cfg = SNAPSHOT_CONFIG[SNAPSHOT_CODE]
    s_init = cfg["s_init"]
    s_end = cfg["s_end"]
    mask_lower = cfg["mask_lower"]
    mask_upper = cfg["mask_upper"]
    interacoes = cfg["interacoes"]
    x_init = cfg["x_init"]
    x_end = cfg["x_end"]
    max_kT = cfg["max_kT"]
    max_dx = cfg["max_dx"]
except KeyError:
    raise ValueError(f"Configuração não encontrada para SNAPSHOT_CODE = {SNAPSHOT_CODE}")

for i in range(len(lines)):
    if (not FULL_PLOT):
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
    if FULL_PLOT:
        lim_sup = 2000
        lim_inf = -2000
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
        length = x_end - x_init
        for j in range(interacoes):
            x1 = x_init + (j * length / interacoes)
            x2 = x_init + ((j + 1) * length / interacoes)
            cond = (x > x1) & (x < x2) & limit_yz

            u = np.mean(u_tot[cond])
            kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV

            x_plot.append((x1 + x2) / 2)
            kT_plot.append(kT)
    kT_plot = fill_nan_nearest(kT_plot)
    kT_plot_filtrado = savgol_filter(kT_plot, window_length=11, polyorder=3)

    if PLOT_INFOS:
        # Medindo a descontinuidade de temperatura
        sep = 18 # numero de pontos entre T1/T2 e a descontinuidade

        x_plot = np.array(x_plot)
        mask_range = (x_plot > mask_lower) & (x_plot < mask_upper) # Considerar só valores dentro desse range para determinar a descontinuidade
        grad_kT = np.full_like(kT_plot_filtrado, 0, dtype=float)
        grad_kT[mask_range] = np.gradient(kT_plot_filtrado[mask_range], x_plot[mask_range])

        # idx = np.where(grad_kT == np.min(grad_kT))[0][0]
        valid_idx = np.where(mask_range)[0]
        if last_snapshot_descontinuity == 0:
            idx = valid_idx[np.argmin(grad_kT[valid_idx])]
        else:
            dx = np.abs(x_plot[valid_idx] - last_snapshot_descontinuity)
            close_idx = valid_idx[dx < max_dx]
            if len(close_idx) > 0:
                idx = close_idx[np.argmin(grad_kT[close_idx])]
            else:
                idx = valid_idx[np.argmin(grad_kT[valid_idx])]
        last_snapshot_descontinuity = x_plot[idx]

        idx_menor = idx + sep
        idx_maior = np.where(kT_plot_filtrado == np.max(kT_plot_filtrado))[0][0]

        T1 = kT_plot_filtrado[idx_menor]
        T2 = kT_plot_filtrado[idx_maior]

        # Calculado o erro de mach pelo filtro Savitzky-Golay
        erro = erro_mach(kT_plot, x_plot, sep, mask_range)
        erro_savgol_mach['time'].append(time)
        erro_savgol_mach['erro'].append(erro)

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
        # plt.close()

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
    ax.plot(x_plot, kT_plot_filtrado, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperatura das partículas')
    # ax.set_title(f"Temperature x Radius - {time:.3f} Gyr")
    ax.set_ylabel('$kT$ (keV)', fontsize=FONT_SIZE)
    ax.set_xlabel('$x$ (kpc)', fontsize=FONT_SIZE)
    ax.tick_params(axis='both', labelsize=FONT_SIZE)
    if FULL_PLOT:
        ax.set_xlim(lim_inf, lim_sup)
        ax.axhline(y=max(kT_plot_filtrado), color='g', linestyle='--', alpha=0.5, 
                   label=f'Maior temperatura: {max(kT_plot_filtrado):.2f}')
    else:
        ax.set_xlim(x_init, x_end)
        ax.set_ylim(0, max_kT)
    if PLOT_INFOS:
        # ax.axhline(y=kT_plot_filtrado[idx_maior], color='g', linestyle='--', alpha=0.2, label=f'T2')
        # ax.axhline(y=kT_plot_filtrado[idx_menor], color='purple', linestyle='--', alpha=0.2, label=f'T1')
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
        os.makedirs(f"{IMAGE_PATH}full_plot", exist_ok=True)
        plt.savefig(f"{IMAGE_PATH}full_plot/t-r_{i:03d}.png")
    else:
        plt.savefig(f"{IMAGE_PATH}t-r_{i:03d}.png")
    plt.close()

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
    plt.close()

    # Print dos outputs
    print(f"Velocidade da onda: {a*0.9778:.2f} km/s")
    print(f"Velocidade do som não chocado: {np.mean(velocities['cs']):.2f} km/s")
    print(f"Velocidade das particulas não chocadas: {np.mean(velocities['u']):.2f} km/s")

    mach_v = (a*0.9778 - np.mean(velocities['u'])) / np.mean(velocities['cs'])
    print(f"Mach pela velocidade: {mach_v:.2f} \n")

    for i in range(len(machs["time"])):
        print(f"Mach por T2/T1 em {machs['time'][i]:.3f} Gyr: {machs['mach'][i]:.2f} ± {erro_savgol_mach['erro'][i]:.2f}")
    print(f"Média dos machs por T2/T1: {np.mean(machs['mach']):.2f} ± {np.mean(erro_savgol_mach['erro']):.2f}")

# ffmpeg -framerate 6 -i plot/t-r_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

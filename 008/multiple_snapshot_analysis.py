import os
import h5py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from functions import *

# CONFIGURAÇÕES POR SNAPSHOT
SNAPSHOT_CONFIG = {
    "0005_0005_0": {
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
    "0005_0005_1000": {
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
    "0005_0005_3000": {
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
    "0005_0006_0": {
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
    "0005_0006_100": {
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
    "0005_0006_1000": {
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
    "0005_0006_2000": {
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
    "0005_0006_3000": {
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
    "0005_0007_0": {
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
    "0005_0007_1000": {
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
    "0005_0007_3000": {
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

FONT_SIZE = 14
IMAGE_NAME = "0005_0006"

SNAPSHOTS = ["0005_0006_0", 
             "0005_0006_100", 
             "0005_0006_1000", 
             "0005_0006_2000", 
             "0005_0006_3000"]

mi = 0.6 # Average molecular weight
Mh = 1.67262192 * 10**(-27) # Proton mass in kg
gamma_var = 5/3 # coeficiente de dilatação adiabática do gás

# cores
pontos = '#4772FF'
linha = '#BACBFF'

plot_data = {"v0": [], "mach": [], "temperature": []}

for SNAPSHOT in SNAPSHOTS:
    SNAPSHOT_CODE = SNAPSHOT
    SNAPSHOT_PATH = "/mnt/d/UFPR/IC/snapshots/" + SNAPSHOT_CODE

    dxdt = {"x": [], "t": []} # objeto para plotar a posição da descontinuidade pelo tempo
    velocities = {"cs": [], "u": []} # objeto para calcular o número de mach pela velocidade
    machs = {"time": [], "mach": []} # objeto para guardar os machs calculados por T2/T1
    erro_savgol_mach = {"time": [], "erro": []} # objeto para guardar os erros dos machs pelo filtro Savitzky-Golay
    last_snapshot_descontinuity = 0 # pega a posição da descontinuidade no ultimo snapshot
    temperatures = []
    lines = open(SNAPSHOT_PATH + "/snapshot.txt", "r").readlines()

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
        if (i < s_init) | (i > s_end):
            continue

        print(f"{bcolors.OKGREEN}{SNAPSHOT} -> Snapshot {i:03d} of {len(lines) - 1:03d}...{bcolors.ENDC}")
        snapshot = f"{SNAPSHOT_PATH}/snapshot_{i:03d}.hdf5"
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
        temperatures.append(T2)

        # Calculado o erro de mach pelo filtro Savitzky-Golay
        erro = erro_mach(kT_plot, x_plot, sep, mask_range)
        erro_savgol_mach['time'].append(time)
        erro_savgol_mach['erro'].append(erro)

        # Calculando as velocidades no gás não chocado
        plot_limit = limit_yz & (x > 200)
        descontinuity_limit = limit_yz & (x > x_plot[idx])
        x_nao_chocado = x_plot[idx] + 250
        volume_limit = limit_yz & (x > x_nao_chocado -100) & (x < x_nao_chocado + 100)

        u = np.mean(u_tot[volume_limit]) # Energia interna no gás não chocado
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV
        cs = velocidade_som_keV(kT) # velocidade do som no gás não chocado
        u = np.mean(vx[volume_limit]) # velocidade do gás não chocado

        velocities['cs'].append(cs)
        velocities['u'].append(u)

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

    # Calcula o Mach pela velocidade
    a, b = np.polyfit(dxdt['t'], dxdt['x'], 1) 
    mach_v = (a*0.9778 - np.mean(velocities['u'])) / np.mean(velocities['cs'])

    plot_data["v0"].append(SNAPSHOT.split("_")[2])
    plot_data["mach"].append(mach_v)
    plot_data["temperature"].append(np.mean(temperatures))

# PLOT
# Magens das escalas
temp = np.array(plot_data["temperature"], dtype=float)
mach = np.array(plot_data["mach"], dtype=float)

temp_margin = 0.05 * (temp.max() - temp.min())
mach_margin = 0.05 * (mach.max() - mach.min())

temp_min = temp.min() - temp_margin
temp_max = temp.max() + temp_margin
mach_min = mach.min() - mach_margin
mach_max = mach.max() + mach_margin

# Plot 
fig, ax1 = plt.subplots(figsize=(8, 6))

# --- Temperatura (eixo Y esquerdo) ---
ax1.plot(
    plot_data["v0"],
    plot_data["temperature"],
    'o-',
    color=pontos,
    label='Temperatura'
)
ax1.set_xlabel('$v_0$ (km/s)', fontsize=FONT_SIZE)
ax1.set_ylabel('$kT$ (keV)', fontsize=FONT_SIZE, color=pontos)
ax1.tick_params(axis='y', labelcolor=pontos)
ax1.set_ylim(temp_min, temp_max)

# --- Mach (eixo Y direito) ---
ax2 = ax1.twinx()
ax2.plot(
    plot_data["v0"],
    plot_data["mach"],
    's--',
    color='tab:red',
    label='Mach'
)
ax2.set_ylabel('Mach', fontsize=FONT_SIZE, color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_ylim(mach_min, mach_max)

# --- Legenda combinada ---
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=FONT_SIZE)

ax1.tick_params(axis='x', labelsize=FONT_SIZE)

plt.tight_layout()
plt.savefig(f"plots/temperature_mach_v0_{IMAGE_NAME}.png")
plt.close()

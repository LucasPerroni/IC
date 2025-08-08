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
gamma_var = 5/3 # coeficiente de dilatação adiabática do gás

# cores
pontos = '#4772FF'
linha = '#BACBFF'

# parâmetros para plot de 4 imagens
# plot_index = 0
# fig, axs = plt.subplots(2, 2, figsize=(16, 12))

dxdt = {"x": [], "t": []}
velocities = {"cs": [], "u": []}
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    if i < 8 or i > 11:
    # if i != 10:
        continue

    print(f"{bcolors.OKGREEN}Plotting snapshot {i:03d} of {len(lines) - 1:03d}...{bcolors.ENDC}")
    snapshot = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    s = h5py.File(snapshot, 'r')

    time = s['Header'].attrs[u'Time'] # Gyr
    mass = s['PartType0']['Masses'][:] # 1e10 Msol
    u_tot = s['PartType0']['InternalEnergy'][:] * 10**6 # (km/s)² -> (m/s)²

    vx = s['PartType0']['Velocities'][:,0]
    vy = s['PartType0']['Velocities'][:,1]
    vz = s['PartType0']['Velocities'][:,2]
    V = np.sqrt(vx**2 + vy**2 + vz**2)

    x = s['PartType0']['Coordinates'][:,0] # kpc
    y = s['PartType0']['Coordinates'][:,1] # kpc
    z = s['PartType0']['Coordinates'][:,2] # kpc
    R = np.sqrt(x**2 + y**2 + z**2) # radius

    x_plot = []
    kT_plot = []
    limit_yz = (y > -100) & (y < 100) & (z > -100) & (z < 100)

    # interacoes = 200
    # length = max(x) - min(x)
    # for j in range(interacoes):
    #     x1 = min(x) + (j * length / interacoes)
    #     x2 = min(x) + ((j + 1) * length / interacoes)
    #     cond = np.argwhere((x > x1) & (x < x2) & limit_yz)

    #     u = np.sum(u_tot[cond]) / len(u_tot[cond])
    #     kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15

    #     x_plot.append((x1 + x2) / 2)
    #     kT_plot.append(kT)

    interacoes = 200
    length = 1000 - 300
    for j in range(interacoes):
        x1 = 300 + (j * length / interacoes)
        x2 = 300 + ((j + 1) * length / interacoes)
        cond = np.argwhere((x > x1) & (x < x2) & limit_yz)

        u = np.sum(u_tot[cond]) / len(u_tot[cond])
        kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV

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

    # Calculando as velocidades no gás não chocado
    plot_limit = limit_yz & (x > 0)
    descontinuity_limit = limit_yz & (x > x_plot[idx])
    x_nao_chocado = sorted(x[descontinuity_limit])[int(len(x[descontinuity_limit])/1.05)]
    volume_limit = limit_yz & (x > x_nao_chocado -100) & (x < x_nao_chocado + 100)

    u = np.mean(u_tot[volume_limit]) # Energia interna no gás não chocado
    kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV
    cs = velocidade_som_keV(kT) # velocidade do som no gás não chocado
    u = np.mean(vx[volume_limit]) # velocidade do gás não chocado

    velocities['cs'].append(cs)
    velocities['u'].append(u)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x[plot_limit], vx[plot_limit], ".")
    ax.set_title(f"Velocity x Radius - {time:.3f} Gyr")
    ax.set_ylabel("velocity (km/s)")
    ax.set_xlabel("x (kpc)")
    ax.axvline(x=x_nao_chocado, color="y", linestyle="--", alpha=0.5, label="Não Chocado")
    ax.axvline(x=x_plot[idx], color="black", linestyle="--", alpha=0.5, label="Descontinuidade")
    ax.axhline(y=-1000, color='r', linestyle='--', alpha=0.5, label="-1000 km/s")
    ax.set_xlim(0, 3000)
    ax.set_ylim(-3000, 3000)
    ax.legend()
    plt.savefig(f"plot/x_vx_{i:03d}.png")

    # Calculando número de Mach
    M = symbols('M')
    eq = Eq((5*M**4 + 14*M**2 - 3) / (16*M**2), T2/T1)
    mach = solve(eq)
    mach = max([m for m in mach if im(m) == 0])
    # print(f"{bcolors.OKCYAN}Mach T2/T1:{bcolors.ENDC} {mach:.4f}")

    dxdt["x"].append(x_plot[idx])
    dxdt['t'].append(time)

    # Plot individual
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    ax.set_title(f"Temperature x Radius - {time:.3f} Gyr - Mach: {mach:.4f}")
    ax.set_ylabel('$kT$ ($keV$)')
    ax.set_xlabel('$x$ ($kpc$)')
    # ax.set_xlim(-3000, 3000)
    ax.set_xlim(300, 900)
    ax.set_ylim(0, 40)
    ax.set_aspect('auto')
    ax.axhline(y=kT_plot[idx_maior], color='g', linestyle='--', alpha=0.5, label=f'T2')
    ax.axhline(y=kT_plot[idx_menor], color='purple', linestyle='--', alpha=0.5, label=f'T1')
    ax.axvline(x=x_plot[idx], color='r', linestyle='--', alpha=0.5, label=f'Descontinuidade')
    ax.legend()
    plt.savefig(f"plot/t-r_{i:03d}.png")

    # Plot para 4 imagens
    # i_plot = plot_index // 2
    # j_plot = plot_index % 2
    # axs[i_plot][j_plot].plot(x_plot, kT_plot, '.', ms = 4, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
    # axs[i_plot][j_plot].set_title(f"Temperature x Radius - {time:.3f} Gyr - Mach: {mach:.4f}")
    # axs[i_plot][j_plot].set_ylabel('$kT$ ($keV$)')
    # axs[i_plot][j_plot].set_xlabel('$x$ ($kpc$)')
    # axs[i_plot][j_plot].set_xlim(300, 900)
    # axs[i_plot][j_plot].set_ylim(0, 40)
    # axs[i_plot][j_plot].set_aspect('auto')
    # axs[i_plot][j_plot].axvline(x=x_plot[idx], color='r', linestyle='--', alpha=0.5, label=f'Descontinuidade')
    # # axs[i_plot][j_plot].axvline(x=x_nao_chocado, color='black', linestyle='--', alpha=0.5, label=f'Não Chocado')
    # axs[i_plot][j_plot].axhline(y=kT_plot[idx_maior], color='g', linestyle='--', alpha=0.5, label=f'T2')
    # axs[i_plot][j_plot].axhline(y=kT_plot[idx_menor], color='purple', linestyle='--', alpha=0.5, label=f'T1')
    # # axs[i_plot][j_plot].axvline(x=x_plot[idx_maior], color='g', linestyle='--', alpha=0.5, label=f'T2')
    # # axs[i_plot][j_plot].axvline(x=x_plot[idx_menor], color='purple', linestyle='--', alpha=0.5, label=f'T1')
    # axs[i_plot][j_plot].legend()
    # plot_index += 1

# save figure - para plot de 4 imagens
# plt.savefig(f"plot/t-r_all.png")

# Ajuste linear a curva de posição pelo tempo
a, b = np.polyfit(dxdt['t'], dxdt['x'], 1)
x_ajuste = np.linspace(min(dxdt['t']), max(dxdt['t']), 100)
y_ajuste = a*x_ajuste + b

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(dxdt['t'], dxdt["x"], ".", label="Posição da descontinuidade")
ax.plot(x_ajuste, y_ajuste, color="lightblue", label="Ajuste linear")
ax.set_xlabel("t (Gyr)")
ax.set_ylabel("x (kpc)")
ax.set_title(f"Posição da Descontinuidade x Tempo - Velocidade: {a:.2f} km/s")
plt.savefig("plot/pos-time.png")

# Print dos outputs
print(f"Velocidade da onda: {a:.2f} km/s")
print(f"Velocidade do som não chocado: {np.mean(velocities['cs']):.2f} km/s")
print(f"Velocidade das particulas não chocadas: {np.mean(velocities['u']):.2f} km/s")

mach_v = (a - np.mean(velocities['u'])) / np.mean(velocities['cs'])
print(f"Mach pela velocidade: {mach_v:.2f}")


# ffmpeg -framerate 6 -i plot/t-r_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

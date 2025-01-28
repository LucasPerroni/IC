import numpy as np
import matplotlib.pyplot as plt

# model = "galaxy"
model = "cluster"
Gamma_h = 1
Gamma_g = 0

# read txt file
p, r = np.loadtxt('txt/density_radius.txt', unpack=True)
p_c, r_c = np.loadtxt('txt/density_radius_comparison.txt', unpack=True)

if model == "cluster":
  p_log, r_log = np.loadtxt('txt/density_radius_log.txt', unpack=True)
  p_g_log, r_g_log = np.loadtxt('txt/density_radius_g_log.txt', unpack=True)
  p_g_c, r_g_c = np.loadtxt('txt/density_radius_g_comparison.txt', unpack=True)

# cores
pontos = '#4772FF'
linha = '#BACBFF'

# legenda
if model == "galaxy":
  funcao = "Exponential Density Function"
elif model == "cluster":
  if Gamma_h == 0:
    funcao = "Dehnen density profile"
  elif Gamma_h == 1:
    funcao = "Hernquist density profile"

  if Gamma_g == 0:
    funcao_g = "Dehnen density profile"
  elif Gamma_g == 1:
    funcao_g = "Hernquist density profile"

if model == "galaxy":
  fig, ax1 = plt.subplots(figsize=(8, 6))
elif model == "cluster":
  fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 6))

  # ax2 - gas density x radius
  ax2.plot(r_g_log, p_g_log, '.', ms = 12, mec = pontos, mfc = pontos, label='Density from Snapshot')
  ax2.plot(r_g_c, p_g_c, '-', ms = 1.5, color = linha, label=funcao_g)
  ax2.set_aspect('auto')
  ax2.set_yscale("log")
  ax2.set_xscale("log")
  ax2.set_xlabel('$r$ ($kpc$)')
  ax2.set_ylabel('$\\rho$ ($M_{\odot} kpc^{-3}$)')
  ax2.set_title('Density x Radius - Gas')
  ax2.legend()

# ax1 - disk/halo density x radius
if model == "galaxy":
  ax1.plot(r, p, '.', ms = 12, mec = pontos, mfc = pontos, label='Density from Snapshot')
  ax1.set_title('Density x Radius - Disk')
  ax1.set_ylabel('$\\rho$ ($M_{\odot} kpc^{-2}$)')
elif model == "cluster":
  ax1.plot(r_log, p_log, '.', ms = 12, mec = pontos, mfc = pontos, label='Density from Snapshot')
  ax1.set_xscale("log")
  ax1.set_title('Density x Radius - Halo')
  ax1.set_ylabel('$\\rho$ ($M_{\odot} kpc^{-3}$)')

ax1.plot(r_c, p_c, '-', ms = 1.5, color = linha, label=funcao)
ax1.set_aspect('auto')
ax1.set_yscale("log")
ax1.set_xlabel('$r$ ($kpc$)')
ax1.legend()

# save figure
plt.savefig('plot/d-r_tests.pdf')

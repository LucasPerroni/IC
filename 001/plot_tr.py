import numpy as np
import matplotlib.pyplot as plt

# cores
pontos = '#4772FF'
linha = '#BACBFF'

# dados e tamanho da figura
keV, r = np.loadtxt('txt/temperature_radius.txt', unpack=True)

fig, ax1 = plt.subplots(figsize=(8, 6))

# ax1 - temperature x radius
ax1.plot(r, keV, '.', ms = 12, mec = pontos, mfc = pontos, label='Temperature from Snapshot')
ax1.set_title('Temperature x Radius')
ax1.set_ylabel('$kT$ ($keV$)')
ax1.set_xlabel('$r$ ($kpc$)')
ax1.set_aspect('auto')
ax1.set_ylim(0, 14)
ax1.legend()

# save figure
plt.savefig('plot/t-r_0000.pdf')

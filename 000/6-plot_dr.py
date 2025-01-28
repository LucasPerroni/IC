import numpy as np
import matplotlib.pyplot as plt

# read txt file
p, r = np.loadtxt('density-radius.txt', unpack=True)

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

# plot data
ax1.plot(r, p, '.')

# aspect ratio
ax1.set_aspect('auto')

# axis ranges
ax1.set_xlim(0, 30)
#ax1.set_ylim(0, 210000000)

ax1.set_yscale("log")

# axis labels
ax1.set_xlabel('$r$ (kpc)')
ax1.set_ylabel('$p$ (Msol/kpc^3)')

# save figure
plt.savefig('plot-dr.pdf')

import numpy as np
import matplotlib.pyplot as plt

# read txt file
x, y, z = np.loadtxt('disk.txt', unpack=True)

# subsample
step = 5
x = x[::step]
y = y[::step]
z = z[::step]

# create figure
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(7, 3.5))

# plot data
ax1.plot(x, y, '.', markersize=0.5)
ax2.plot(x, z, '.', markersize=0.5)

# aspect ratio
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# axis ranges
ax1.set_xlim(-20, 20)
ax1.set_ylim(-20, 20)
ax2.set_xlim(-20, 20)
ax2.set_ylim(-20, 20)

# axis labels
ax1.set_xlabel('$x$ (kpc)')
ax1.set_ylabel('$y$ (kpc)')
ax2.set_xlabel('$x$ (kpc)')
ax2.set_ylabel('$z$ (kpc)')

# adjust space between subplots
fig.subplots_adjust(hspace=0, wspace=0.5)

# save figure
plt.savefig('plot.pdf')

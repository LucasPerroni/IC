import numpy as np
import h5py
import math

# path to the snapshot
snapshot = '../ic/snapshot_0000.hdf5'

# load snapshot
s = h5py.File(snapshot, 'r')

# read time
time = s['Header'].attrs[u'Time']

# read masses of the disk particles
mass = s['PartType2']['Masses'][:]

# read coordinates of the disk particles
x = s['PartType2']['Coordinates'][:,0]
y = s['PartType2']['Coordinates'][:,1]
z = s['PartType2']['Coordinates'][:,2]

# cylindrical radius
R = np.sqrt(x**2 + y**2 + z**2)

interacoes = 20
rmax = np.max(R)
passo = rmax/interacoes
r1 = 0
r2 = passo

densidades = []
raios = []

for i in range(0, interacoes):
	cond = np.argwhere(np.logical_and(R > r1, R < r2))
	m = np.sum(mass[cond])
	vol = (4/3 * math.pi * r2**3) - (4/3 * math.pi * r1**3)
	densidade = m * 1e10 / vol
	raio = (r2 + r1) / 2
	densidades.append(densidade)
	raios.append(raio)
	r1 = r1 + passo
	r2 = r2 + passo

# open file for writing
f = open('density-radius.txt', 'w')

# write file header (column names and units)
f.write('# p(Msol/kpc^3) y(kpc)\n')

# write data into file
for i in range(0, interacoes):
    f.write('%f %f \n' % (densidades[i], raios[i]))

f.close()

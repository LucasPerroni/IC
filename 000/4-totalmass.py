import numpy as np
import h5py

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

#total disk mass
Mtot = np.sum(mass)

#units of Msun
Mtot = Mtot * 1e10

print('Mtot = %.2e Msun' % Mtot)

# cylindrical radius
R = np.sqrt(x**2 + y**2)

#indices of particles within 6 kpc
cond = np.argwhere( R < 6.0 )

#sum of masses of particles that satisfy the condition
Min = np.sum( mass[cond] )

#units of Msun
Min = Min * 1e10

print('M(R<6kpc) = %.2e Msun' % Min)

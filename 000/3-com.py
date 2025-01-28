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

# compute centre of mass
xcom = np.sum(x*mass) / np.sum(mass)
ycom = np.sum(y*mass) / np.sum(mass)
zcom = np.sum(z*mass) / np.sum(mass)

print('time = %f Gyr' % time)
print('xcom = %f kpc' % xcom)
print('ycom = %f kpc' % ycom)
print('zcom = %f kpc' % zcom)

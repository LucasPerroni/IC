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


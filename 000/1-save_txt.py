import h5py

# path to the snapshot
snapshot = '../ic/snapshot_0000.hdf5'

# load snapshot
s = h5py.File(snapshot, 'r')

# read coordinates of the disk particles
x = s['PartType2']['Coordinates'][:,0]
y = s['PartType2']['Coordinates'][:,1]
z = s['PartType2']['Coordinates'][:,2]

# length of array
Ndisk = len(x)

# open file for writing
f = open('disk.txt', 'w')

# write file header (column names and units)
f.write('# x(kpc) y(kpc) z(kpc)\n')

# write data into file
for i in range(0, Ndisk):
    f.write('%f %f %f \n' % (x[i], y[i], z[i]))

f.close()

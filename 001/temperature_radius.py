import math
import h5py
import numpy as np

snapshot = "../../hdf5/snapshot_0000.hdf5"
s = h5py.File(snapshot, 'r')

# read time
time = s['Header'].attrs[u'Time']

mass = s['PartType0']['Masses'][:] # 1e10 Msol
u_tot = s['PartType0']['InternalEnergy'][:] * 10**6 # (km/s)² -> (m/s)²

x = s['PartType0']['Coordinates'][:,0] # kpc
y = s['PartType0']['Coordinates'][:,1] # kpc
z = s['PartType0']['Coordinates'][:,2] # kpc

mi = 0.6 # Average molecular weight
Mh = 1.67262192 * 10**(-27) # Proton mass in kg
R = np.sqrt(x**2 + y**2 + z**2) # radius

# number of points in the plot
interacoes = 20

# maximum radius
rmax = np.max(R)

# low and high radius for temperature calc, as well as the function step size
passo = rmax/interacoes
r1 = 0
r2 = passo

# arrays for temperature and radius points
temperaturas = []
raios = []

for i in range(interacoes):
  # get radius between r1 and r2 and the internal energies in those points
  cond = np.argwhere(np.logical_and(R > r1, R < r2))
  u = np.sum(u_tot[cond]) / len(u_tot[cond]) # Average energy density

  kT = (u * (2 * mi * Mh) / 3) * 6.241506 * 10**15 # Temperature in J -> keV
  raio = (r2 + r1) / 2

  # append values to arrays
  temperaturas.append(kT)
  raios.append(raio)

  # update radius
  r1 += passo
  r2 += passo


# save data in txts
name = 'txt/temperature_radius.txt' # linear values for temperature and radius
size = len(temperaturas)
f = open(name, 'w')

f.write('# kT(keV) r(kpc) \n')
for i in range(0, size):
  f.write('%f %f \n' % (temperaturas[i], raios[i]))
f.close()

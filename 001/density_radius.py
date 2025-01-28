import math
import h5py
import numpy as np

snapshot = '../../hdf5/snapshot_0002.hdf5'
s = h5py.File(snapshot, 'r')

# model = "galaxy"
model = "cluster"

# Galaxy disk variables
Md = 5 # stellar disk mass
Rd = 3.5 # stellar disk radial scale factor

# Cluster halo variables
Mh = 50000 # halo mass
Ah = 200 # scale factor in the Dehnen density profile
Gamma_h = 1 # Dehnen density profile free parameter. If gamma = 1, then the profile
            # is equal to a Hernquist profile, featuring a central cusp. If gamma = 0,
            # then the profile features a central core.

# Cluster gas variables
Mg = 15000 # gas mass
Ag = 200 # scale factor in the Dehnen density profile
Gamma_g = 0

# GALAXIES:
#   PartType0: Gas
#   PartType1: Halo
#   PartType2: Disk
#   PartType3: Bulge

# CLUSTERS:
#   PartType0:  Gas
#   PartType1:  Dark Matter Halo

# read time
time = s['Header'].attrs[u'Time']

# read masses of the disk particles
if model == "galaxy":
  mass = s['PartType2']['Masses'][:] # Galaxy
elif model == "cluster":
  mass = s['PartType1']['Masses'][:] # Halo
  mass_g = s['PartType0']['Masses'][:] # Gas

# read coordinates of the disk particles
if model == "galaxy":
  x = s['PartType2']['Coordinates'][:,0] # Galaxy
  y = s['PartType2']['Coordinates'][:,1] # Galaxy
  z = s['PartType2']['Coordinates'][:,2] # Galaxy
elif model == "cluster":
  x = s['PartType1']['Coordinates'][:,0] # Halo
  y = s['PartType1']['Coordinates'][:,1] # Halo
  z = s['PartType1']['Coordinates'][:,2] # Halo

  x_g = s['PartType0']['Coordinates'][:,0] # Gas
  y_g = s['PartType0']['Coordinates'][:,1] # Gas
  z_g = s['PartType0']['Coordinates'][:,2] # Gas

# radius
if model == "galaxy":
  R = np.sqrt(x**2 + y**2) # for galaxies
elif model == "cluster":
  R = np.sqrt(x**2 + y**2 + z**2) # halo
  R_g = np.sqrt(x_g**2 + y_g**2 + z_g**2) # gas

# number of points in the plot
interacoes = 20

# maximum radius
rmax = np.max(R)
rmax_g = np.max(R_g)

# low and high radius for density calc, as well as the function step size
passo = rmax/interacoes
r1 = 0
r2 = passo

# arrays for density and radius points
densidades = []
raios = []

for i in range(interacoes):
  # get radius between r1 and r2 and the masses in those points
  cond = np.argwhere(np.logical_and(R > r1, R < r2))
  m = np.sum(mass[cond])

  # calculate density and radius
  if model == "galaxy":
    vol = (math.pi * r2**2) - (math.pi * r1**2) # for galaxies
  elif model == "cluster":
    vol = (4/3 * math.pi * r2**3) - (4/3 * math.pi * r1**3) # for clusters

  densidade = m * 1e10 / vol
  raio = (r2 + r1) / 2

  # append values to arrays
  densidades.append(densidade)
  raios.append(raio)

  # update radius
  r1 = r1 + passo
  r2 = r2 + passo

if model == "cluster":
  # Define the number of bins for the logarithmic radial bins
  n_bins = 21

  # logarithmically spaced radial bins
  log_bins = np.logspace(np.log10(1), np.log10(rmax), n_bins)
  log_g_bins = np.logspace(np.log10(1), np.log10(rmax_g), n_bins)

  # arrays for density and radius points
  densidades_log = []
  densidades_g_log = []
  raios_log = []
  raios_g_log = []

  for i in range(n_bins - 1):
    # low and high radius for density calc in logarithmically spaced radial bins
    r1_log = log_bins[i]
    r2_log = log_bins[i + 1]

    r1_g_log = log_g_bins[i]
    r2_g_log = log_g_bins[i + 1]

    # get radius between r1 and r2 and the masses in those points
    cond = np.argwhere(np.logical_and(R > r1_log, R < r2_log))
    m = np.sum(mass[cond])

    cond_g = np.argwhere(np.logical_and(R_g > r1_g_log, R_g < r2_g_log))
    m_g = np.sum(mass_g[cond_g])

    # calculate density and radius
    vol = (4/3 * math.pi * r2_log**3) - (4/3 * math.pi * r1_log**3) # for clusters
    vol_g = (4/3 * math.pi * r2_g_log**3) - (4/3 * math.pi * r1_g_log**3) # for clusters

    densidade = m * 1e10 / vol
    raio = (r2_log + r1_log) / 2

    densidade_g = m_g * 1e10 / vol_g
    raio_g = (r2_g_log + r1_g_log) / 2

    # append values to arrays
    densidades_log.append(densidade)
    raios_log.append(raio)

    densidades_g_log.append(densidade_g)
    raios_g_log.append(raio_g)

# Analytic function data
min_radius = 0.1 # first point of the array
step = 0.1

densidades_comp = []
densidades_g_comp = []
raios_comp = []
raios_g_comp = []

for r in np.arange(min_radius, rmax, step):
  if model == "galaxy":
    p = Md * 1e10 / (2 * math.pi * Rd**2) * math.e**(-r / Rd)
  elif model == "cluster":
    p = ((3 - Gamma_h) * Mh * 1e10 / (4 * math.pi)) * (Ah / (r**Gamma_h * (r + Ah)**(4 - Gamma_h)))

  densidades_comp.append(p)
  raios_comp.append(r)

if model == "cluster":
  for r_g in np.arange(min_radius, rmax_g, step):
    p_g = ((3 - Gamma_g) * Mg * 1e10 / (4 * math.pi)) * (Ag / (r_g**Gamma_g * (r_g + Ag)**(4 - Gamma_g)))

    densidades_g_comp.append(p_g)
    raios_g_comp.append(r_g)


# save data in txts
name = 'txt/density_radius.txt' # linear values for density and radius
name_c = 'txt/density_radius_comparison.txt' # values from analytic function
name_log = 'txt/density_radius_log.txt' # log values for density and radius
name_g_log = 'txt/density_radius_g_log.txt' # log values for density and radius of the gas
name_g_c = 'txt/density_radius_g_comparison.txt' # values from analytic function of the gas

size = len(densidades)
size_comp = len(densidades_comp)

# open file for writing
f = open(name, 'w')
f_c = open(name_c, 'w')

# write file header (column names and units)
if model == "galaxy":
  f.write('# p(Msol/kpc^2) r(kpc) \n')
  f_c.write('# p(Msol/kpc^2) r(kpc) \n')
elif model == "cluster":
  f.write('# p(Msol/kpc^3) r(kpc) \n')
  f_c.write('# p(Msol/kpc^3) r(kpc) \n')

# write data into file
for i in range(0, size):
  f.write('%f %f \n' % (densidades[i], raios[i]))
f.close()

for i in range(0, size_comp):
  f_c.write('%f %f \n' % (densidades_comp[i], raios_comp[i]))
f_c.close()

if model == "cluster":
  size_log = len(densidades_log)
  size_g_comp = len(densidades_g_comp)

  f_log = open(name_log, 'w')
  f_g_log = open(name_g_log, 'w')
  f_g_c = open(name_g_c, 'w')

  f_log.write('# p(Msol/kpc^3) r(kpc) \n')
  f_g_log.write('# p(Msol/kpc^3) r(kpc) \n')
  f_g_c.write('# p(Msol/kpc^3) r(kpc) \n')


  for i in range(0, size_log):
    f_log.write('%f %f \n' % (densidades_log[i], raios_log[i]))
    f_g_log.write('%f %f \n' % (densidades_g_log[i], raios_g_log[i]))
  f_log.close()
  f_g_log.close()

  for i in range(0, size_g_comp):
    f_g_c.write('%f %f \n' % (densidades_g_comp[i], raios_g_comp[i]))
  f_g_c.close()

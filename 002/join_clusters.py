import h5py
from sys import argv
import numpy as np
import argparse

'''
python join_clusters.py clusterA.hdf5 clusterB.hdf5 -o both.hdf5 [-rP X Y Z] [-rV vX vY vZ]

'''

def main():
    parser = argparse.ArgumentParser(description="Processa arquivos HDF5 com opções configuráveis.")

    # Argumentos posicionais obrigatórios
    parser.add_argument("input_file1", type=str, help="Primeiro arquivo HDF5 de entrada")
    parser.add_argument("input_file2", type=str, help="Segundo arquivo HDF5 de entrada")

    # Argumento opcional -o
    parser.add_argument("-o", "--output", type=str, required=True, help="Arquivo HDF5 de saída")

    # Argumento opcional -rP
    parser.add_argument("-rP", "--relative_position", nargs=3, type=int, help="Parâmetros de intervalo (três inteiros)", metavar=('X', 'Y', 'Z'))
    
    # Argumento opcional -rV
    parser.add_argument("-rV", "--relative_velocity", nargs=3, type=int, help="Parâmetros de intervalo (três inteiros)", metavar=('vX', 'vY', 'vZ'))

    # Parse dos argumentos
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = main()

snapshot_1 = args.input_file1
snapshot_2 = args.input_file2
snapshot_out = args.output

if args.relative_position:
	xshift = args.relative_position[0]
	yshift = args.relative_position[1]
	zshift = args.relative_position[2]
else:
	xshift = 2000
	yshift = 0
	zshift = 0

if args.relative_velocity:
	vxshift = args.relative_velocity[0]
	vyshift = args.relative_velocity[1]
	vzshift = args.relative_velocity[2]
else:
	vxshift = -2000
	vyshift = 0
	vzshift = 0

shift = np.array([xshift, yshift, zshift])
vshift = np.array([vxshift, vyshift, vzshift])

#read snapshot_1
with h5py.File(snapshot_1, 'r') as s1:
    gas_mass1 = s1['PartType0']['Masses'][:]
    gas_pos1 = s1['PartType0']['Coordinates'][:]
    gas_vel1 = s1['PartType0']['Velocities'][:]
    gas_rho1 = s1['PartType0']['Density'][:]
    gas_u1 = s1['PartType0']['InternalEnergy'][:]

    halo_mass1  = s1['PartType1']['Masses'][:]
    halo_pos1 = s1['PartType1']['Coordinates'][:]
    halo_vel1 = s1['PartType1']['Velocities'][:]
    halo_ids1 = s1['PartType1']['ParticleIDs'][:]

#read snapshot_2
with h5py.File(snapshot_2, 'r') as s2:
    gas_mass2 = s2['PartType0']['Masses'][:]
    gas_pos2 = s2['PartType0']['Coordinates'][:]
    gas_vel2 = s2['PartType0']['Velocities'][:]
    gas_rho2 = s2['PartType0']['Density'][:]
    gas_u2 = s2['PartType0']['InternalEnergy'][:]

    halo_mass2  = s2['PartType1']['Masses'][:]
    halo_pos2 = s2['PartType1']['Coordinates'][:]
    halo_vel2 = s2['PartType1']['Velocities'][:]
    halo_ids2 = s2['PartType1']['ParticleIDs'][:]

#concatenate halo (shifting snapshot2)
halo_mass = np.concatenate([halo_mass1, halo_mass2])
halo_pos = np.concatenate([halo_pos1, halo_pos2 + shift])
halo_vel = np.concatenate([halo_vel1, halo_vel2 + vshift])

#concatenate gas (shifting snapshot2)
gas_mass = np.concatenate([gas_mass1, gas_mass2])
gas_pos = np.concatenate([gas_pos1, gas_pos2 + shift])
gas_vel = np.concatenate([gas_vel1, gas_vel2 + vshift])
gas_rho = np.concatenate([gas_rho1, gas_rho2])
gas_u = np.concatenate([gas_u1, gas_u2])

Ngas = len(gas_mass)
Nhalo = len(halo_mass)

Npart = [Ngas, Nhalo, 0, 0, 0, 0]

gas_ids = np.arange(1, Ngas+1, dtype=int)
halo_ids = np.arange(Ngas+1, Ngas+Nhalo+1, dtype=int)
#print(gas_ids)
#print(halo_ids)

#create hdf5 file for writing
with h5py.File(snapshot_out, 'w') as f:
    header = f.create_group('Header')
    header.attrs['NumPart_ThisFile'] = np.asarray(Npart)
    header.attrs['NumPart_Total'] = np.asarray(Npart)
    header.attrs['NumPart_Total_HighWord'] = 0 * np.asarray(Npart)
    header.attrs['MassTable'] = np.zeros(6)
    header.attrs['Time'] = float(0.0)
    header.attrs['Redshift'] = float(0.0)
    header.attrs['BoxSize'] = float(0.0)
    header.attrs['NumFilesPerSnapshot'] = int(1)
    header.attrs['Omega0'] = float(0.0)
    header.attrs['OmegaLambda'] = float(0.0)
    header.attrs['HubbleParam'] = float(1.0)
    header.attrs['Flag_Sfr'] = int(0.0)
    header.attrs['Flag_Cooling'] = int(0)
    header.attrs['Flag_StellarAge'] = int(0)
    header.attrs['Flag_Metals'] = int(0)
    header.attrs['Flag_Feedback'] = int(0)
    header.attrs['Flag_DoublePrecision'] = 0
    header.attrs['Flag_IC_Info'] = 0

    gas = f.create_group(f'PartType0')
    gas.create_dataset('Masses', data=gas_mass, dtype='float32')
    gas.create_dataset('ParticleIDs', data=gas_ids, dtype='uint32')
    gas.create_dataset('Coordinates', data=gas_pos, dtype='float32')
    gas.create_dataset('Velocities', data=gas_vel, dtype='float32')
    gas.create_dataset('Density', data=gas_rho, dtype='float32')
    gas.create_dataset('InternalEnergy', data=gas_u, dtype='float32')
    
    halo = f.create_group(f'PartType1')
    halo.create_dataset('Masses', data=halo_mass, dtype='float32')
    halo.create_dataset('ParticleIDs', data=halo_ids, dtype='uint32')
    halo.create_dataset('Coordinates', data=halo_pos, dtype='float32')
    halo.create_dataset('Velocities', data=halo_vel, dtype='float32')


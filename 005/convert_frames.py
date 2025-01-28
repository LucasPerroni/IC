import subprocess

# Caminho para o diretório contendo as imagens
GAS_DIR = './frames/gas/'
HALO_DIR = "./frames/halo/"
OUTPUT_DIR = "./frames/join/"

PATH_TO_HDF5 = "./hdf5/"
PATH_TO_PLOT = "./plots/"
lines = open(PATH_TO_HDF5 + "snapshot.txt", "r").readlines()

# Itera sobre os pares de imagens e executa o comando 'convert'
for i in range(len(lines)):
    image1 = f"{HALO_DIR}frame_.{i:05d}.jpg" 
    image2 = f"{GAS_DIR}frame_.{i:05d}.jpg"
    output_image = f"{OUTPUT_DIR}out_{i:05d}.jpg"
    
    # Executa o comando 'convert'
    subprocess.run(['convert', image1, image2, '+append', output_image])

print("Conversão concluída!")

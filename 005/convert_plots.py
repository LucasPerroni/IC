import subprocess

# Caminho para o diretório contendo as imagens
FRAME_DIR = './frames/join/'
PLOT_DIR = "./plots/"
OUTPUT_DIR = "./images/"

PATH_TO_HDF5 = "./hdf5/"
PATH_TO_PLOT = "./plots/"
lines = open(PATH_TO_HDF5 + "snapshot.txt", "r").readlines()

# Itera sobre os pares de imagens e executa o comando 'convert'
for i in range(1, len(lines)):
    image1 = f"{FRAME_DIR}out_{i:05d}.jpg"
    image2 = f"{PLOT_DIR}d-r_{i:03d}.jpg" 
    output_image = f"{OUTPUT_DIR}image_{i:05d}.jpg"
    
    # Executa o comando 'convert'
    subprocess.run(['convert', image1, image2, '-append', output_image])

print("Conversão concluída!")

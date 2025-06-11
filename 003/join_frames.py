import os
import subprocess

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Caminho para o diretório contendo as imagens
IMAGE_PATH = "/home/lucasbondep/ic_astronomia/main/003/frames/"
GAS_DIR = f'{IMAGE_PATH}gas/'
DM_DIR = f"{IMAGE_PATH}dm/"
OUTPUT_DIR = f"{IMAGE_PATH}join/"

PATH_TO_HDF5 = "/mnt/d/IC/snapshots/"
# lines = open(PATH_TO_HDF5 + "snapshot.txt", "r").readlines()
# 
# # Itera sobre os pares de imagens e executa o comando 'convert'
# for i in range(len(lines)):
#     print(f"{bcolors.OKGREEN}Joining frames {i}...{bcolors.ENDC}")
    
#     image1 = f"{DM_DIR}dm_{i:03d}.png"
#     image2 = f"{GAS_DIR}gas_{i:03d}.png"
#     output_image = f"{OUTPUT_DIR}join_{i:03d}.png"

#     # Executa o comando 'convert'
#     subprocess.run(['convert', image1, image2, '+append', output_image])

img00 = f'{GAS_DIR}temp_008.png'
img01 = f'{GAS_DIR}temp_009.png'
img10 = f'{GAS_DIR}temp_010.png'
img11 = f'{GAS_DIR}temp_011.png'
output_image = f"{OUTPUT_DIR}join_all.png"

subprocess.run(['convert', '+append', img00, img01, 'linha1.png'])
subprocess.run(['convert', '+append', img10, img11, 'linha2.png'])
subprocess.run(['convert', '-append', 'linha1.png', 'linha2.png', output_image])
os.remove('linha1.png')
os.remove('linha2.png')

print("Conversão concluída!")

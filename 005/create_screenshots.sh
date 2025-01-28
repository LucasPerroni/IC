#!/bin/bash

# Le todos os snapshots e salva seus nomes em um txt
ls hdf5/snapshot* > hdf5/snapshot.txt

# Executa o glnemo2 com os parâmetros desejados para o halo
glnemo2 hdf5/snapshot.txt halo play=t grid=FALSE screenshot=frames/halo/frame_

# Executa o glnemo2 com os parâmetros desejados para o gas
glnemo2 hdf5/snapshot.txt gas play=t grid=FALSE screenshot=frames/gas/frame_

# Executa o script Python
python3 convert_frames.py

# Cria a animacao com os frames
ffmpeg -framerate 8 -i frames/join/out_%05d.jpg -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y animation.mp4

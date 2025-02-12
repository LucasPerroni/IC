#!/bin/bash

rm frames/dm/dm_*
rm frames/gas/gas_*
rm frames/gas/temp_*
rm frames/join/join_*

ls /mnt/d/IC/snapshots/snapshot_* > /mnt/d/IC/snapshots/snapshot.txt

python3 create_frames.py

python3 join_frames.py

ffmpeg -framerate 8 -i frames/join/join_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y dm_gas_rho.mp4
ffmpeg -framerate 8 -i frames/gas/temp_%03d.png -vf "scale=930:748" -c:v libx264 -pix_fmt yuv420p -y gas_temp.mp4

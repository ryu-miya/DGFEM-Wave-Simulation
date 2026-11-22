#!/bin/bash
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96

for ne in 50 100 200 400 800;
do  
    g++ -O2 -DNe=$ne wave.cpp -lm
    ./a.out >> conv.dat
    mkdir data
    mv wave.dat data/wave_$ne.dat
    python snapshot.py $ne
    mkdir movie
    ffmpeg -r 30 -i snapshot/snapshot_$ne/frame_%05d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2"  -vcodec libx264 -profile:v high -preset slower -pix_fmt yuv420p -g 30 -y movie/snapshot_$ne.mp4
done
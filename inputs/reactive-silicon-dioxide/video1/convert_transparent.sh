#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"
set -e

# 1) Generate white background movie
for file in _light.*.ppm; 
do 
	convert $file -transparent white ${file:0:6}.png;
done
img2webp -o ../../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/SiO_gif_light.webp -q 40 -mixed -d 22.22 _light*.png

# 2) Generate black background movie
for file in _dark.*.ppm; 
do 
	convert $file -transparent black ${file:0:5}.png;
done
img2webp -o ../../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/SiO_gif_dark.webp -q 40 -mixed -d 22.22 _dark*.png

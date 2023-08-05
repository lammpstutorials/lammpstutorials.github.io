#!/bin/bash

# assuming _dark.00062.ppm or _light.00062.ppm as file name

export LC_NUMERIC="en_US.UTF-8"
set -e

# 1) Generate white background movie
for file in _light.*.ppm; 
do 
	convert $file -transparent white ${file:0:12}.png;
done
img2webp -o ../../../../../../docs/sphinx/source/tutorials/figures/level1/breaking-a-carbon-nanotube/CNT_gif_light.webp -q 10 -mixed -d 22.22 _light*.png

# 2) Generate black background movie
for file in _dark.*.ppm; 
do 
	convert $file -transparent black ${file:0:11}.png;
done
img2webp -o ../../../../../../docs/sphinx/source/tutorials/figures/level1/breaking-a-carbon-nanotube/CNT_gif_dark.webp -q 10 -mixed -d 22.22 _dark*.png

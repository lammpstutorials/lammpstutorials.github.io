#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"
set -e

# 1) Generate white background movie
for file in _light.*.ppm; 
do 
	echo ${file:0:12}.png
	convert $file -transparent white ${file:0:12}.png;
done
img2webp -o ../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/hydronium_transfert_light.webp -q 40 -mixed -d 22.22 _light*.png

# 2) Generate black background movie
for file in _dark.*.ppm; 
do 
	echo ${file:0:11}.png
	convert $file -transparent black ${file:0:11}.png;
done
img2webp -o ../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/hydronium_transfert_dark.webp -q 40 -mixed -d 22.22 _dark*.png

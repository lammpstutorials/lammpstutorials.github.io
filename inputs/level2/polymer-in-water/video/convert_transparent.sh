#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"
set -e


cd mergePEGH2O/

# 1) Generate white background movie
for file in light.*.ppm; 
do 
	convert $file -crop 500x500+55+15 -transparent white ${file:0:11}.png; # -trim
done

img2webp -o ../../../../../../docs/sphinx/source/tutorials/figures/level2/polymer-in-water/PEG-light.webp -q 30 -mixed -d 50 light*.png

# 2) Generate black background movie
for file in dark.*.ppm; 
do 
	convert $file -crop 500x500+55+15 -transparent black ${file:0:10}.png; # -trim
done

img2webp -o ../../../../../../docs/sphinx/source/tutorials/figures/level2/polymer-in-water/PEG-dark.webp -q 30 -mixed -d 50 dark*.png

cd ..


#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"
set -e

# 1) Generate white background movie
for file in light.*.ppm; 
do 
	convert $file -transparent white ${file:0:11}.png;
done

img2webp -o ../../../../../docs/sphinx/source/tutorials/figures/level3/free-energy-calculation/avatar_light.webp -q 30 -mixed -d 70 light*.png

# 2) Generate black background movie
for file in dark.*.ppm; 
do 
	convert $file -transparent black ${file:0:10}.png;
done

img2webp -o ../../../../../docs/sphinx/source/tutorials/figures/level3/free-energy-calculation/avatar_dark.webp -q 30 -mixed -d 70 dark*.png

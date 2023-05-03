#!/bin/sh

for file in untitled.*.ppm; 
do 
	convert $file -transparent white $file.png;
	#rm $file
done

img2webp -o video-solution-light.webp -q 40 -mixed -d 66.66 *.png

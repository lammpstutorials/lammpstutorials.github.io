#!/bin/sh

for file in untitled.*.ppm; 
do 
	convert $file -transparent black $file.png;
	#rm $file
done

img2webp -o video-solution-dark.webp -q 40 -mixed -d 66.66 *.png

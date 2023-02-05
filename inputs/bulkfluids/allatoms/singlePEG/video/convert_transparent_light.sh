

for file in untitled.*.ppm; 
do 
	convert $file -transparent white $file.png;
	rm $file
done

img2webp -o PEG.webp -q 40 -mixed -d 150 *.png

#ffmpeg -framerate 25 -i untitled.%05d.ppm.png -c copy binary_LJ_fluid.mov
#img2webp -o binary_LJ_fluid.webp -q 40 -mixed -d 66.66 *.png
#ffmpeg -r 60 -i untitled.%05d.ppm.png -vcodec libx264 -crf 0  -pix_fmt yuv420p myvideo.mp4
#ffmpeg -i myvideo.mp4 -vcodec libwebp -filter:v fps=fps=20 -lossless 1 -loop 0 -preset default -an -vsync 0 myvideo.webp
#ffmpeg -i untitled.%05d.ppm.png -vcodec png z.mov
#mpeg -i untitled.mov -vcodec libwebp -filter:v fps=fps=20 -lossless 1 -loop 0 -preset default -an -vsync 0 myvideo.webp

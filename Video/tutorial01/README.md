## How to generate the movie

1. Run the input.lammps using LAMMPS

2. Open the dump.lammpstrj file using VMD. Import the state.vmd file to have the same represenatation.

3. Generate image using the movie maker of VMD (without deleting the file)

4. Eventually reshape the images using

for file in untitled.*.ppm; do convert $file -crop 900x468+0+0 output-$file; done

5. merge the images into a video using

ffmpeg -r 60 -i output-untitled.%05d.ppm -vcodec libx264 -crf 0  -pix_fmt yuv420p myvideo.mp4

6. eventually, reduce the size of the video for web integration

ffmpeg -i myvideo.mp4 -vcodec libwebp -filter:v fps=fps=60 -lossless 1 -loop 0 -preset default -an -vsync 0 myvideo.webp
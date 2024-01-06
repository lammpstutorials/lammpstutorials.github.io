
#!/bin/bash

set -e

for file in untitled.*.ppm; do convert $file -crop 1912x402+0+0 $file; done
ffmpeg -r 30 -i untitled.%05d.ppm -vcodec libx264 -crf 0  -pix_fmt yuv420p myvideo.mp4

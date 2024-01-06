
#!/bin/bash

set -e

ffmpeg -r 30 -i untitled.%05d.ppm -vcodec libx264 -crf 0  -pix_fmt yuv420p myvideo.mp4

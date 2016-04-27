!/bin/bash

ffmpeg -framerate 25 -i frame-%08d.png -c:v libx264 -pix_fmt yuv420p movie.mp4

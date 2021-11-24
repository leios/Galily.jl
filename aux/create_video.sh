#!/bin/bash

ffmpeg -i out%04d.png -vf palettegen palette.png
ffmpeg -r 60 -f image2 -i palette.png -i out%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4

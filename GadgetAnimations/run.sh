N=120
ROT=360

./convert $N
echo "#declare ROT = $ROT;" > rotation.inc
povray ani.pov +KFF$N +KF$N
ffmpeg -y -framerate 30 -i ani%0${#N}d.png -c:v libx264 -preset slow -crf 17 animation.mp4

mkdir frames
mv *png frames/.
mv *mp4 frames/.
rm *.inc



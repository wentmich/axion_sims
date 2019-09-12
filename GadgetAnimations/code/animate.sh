#!/bin/bash
##
#SBATCH --job-name=GadgetAnimation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=6:00:00
#SBATCH --account=bsafdi
#SBATCH --partition=standard
#SBATCH --mail-user=wentmich@umich.edu

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

#  Change to the directory you submitted from
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

 # -----------------------------------------------------------------------

CODE_DIR=/nfs/turbo/bsafdi/wentmich/GadgetAnimations/code-hist
cp $CODE_DIR/ani.pov .
cp $CODE_DIR/convert .
cp $CODE_DIR/povray.ini .

./convert $N
povray ani.pov +KFF$N +KF$N
ffmpeg -r $FPS  -f image2 -s 1920x1080 -i ani%0${#N}d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p animation.mp4

mkdir frames
mv *png frames/.
mv *mp4 frames/.
rm ani.pov convert povray.inv
rm *.inc



#!/bin/bash

GRO=~/mass-storage/data/2014-08-cpTMV/md_7md/mdtraj.gro
TPR=~/mass-storage/data/2014-08-cpTMV/md_7md/md.tpr
RMSD=~/mass-storage/data/2014-08-cpTMV/md_7md/rmsd.xyz
IMG_DIR=$HOME/Dropbox/Public/2014-10-06-activity

#times=(0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0 5.0 7.5 10.0)
##times=(19.9)
#DIR=$HOME/Dropbox/Public/2014-10-06-activity
#mkdir $DIR
#for i in ${times[@]}; do
#    FNAME=angle_act
#    python interface.py -rmsd_grid png -rmsd_dt $i -rmsd_type angle -rmsd_activityparm 1.2 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done
#
#for i in ${times[@]}; do
#    FNAME=rmsd_act
#    python interface.py -rmsd_grid png -rmsd_dt $i -rmsd_type rmsd -rmsd_activityparm .38 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done
#
#for i in ${times[@]}; do
#    FNAME=diff_act
#    python interface.py -rmsd_grid png -rmsd_dt $i -rmsd_type rmsd -rmsd_scaletime -rmsd_activityparm .12 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done


times=(19.9)

for i in ${times[@]}; do
    FNAME=rmsd_act
    python interface.py -rmsd_grid png -rmsd_dt $i -rmsd_type angle -rmsd_activityparm 1.2 20 -rmsd_colorrange .1 .9 -rmsd_grid_file $FNAME
    mv $FNAME.png $IMG_DIR/$FNAME."$i"ps.png
done

#!/bin/bash

GRO=~/mass-storage/data/2014-08-cpTMV/md_7md/mdtraj.gro
TPR=~/mass-storage/data/2014-08-cpTMV/md_7md/md.tpr
RMSD=~/mass-storage/data/2014-08-cpTMV/md_7md/rmsd.xyz
#python RMSDAnalyze.py -rmsd_out $RMSD -rmsd_dt 10
#trjconv -b 10 -f $GRO -s $TPR -o local_traj.gro
#rm activity.xyz
#for dt in {1..10}; do
#    python RMSDAnalyze.py -act_out activity.xyz -act_dt $dt -act_frames 1 -act_dist $(echo "$dt * .1" | bc) -act_slope $(echo "10. / $dt * 3." | bc -l) -act_openflag 'a'
#done
#python RMSDAnalyze.py -gro_copy $GRO grocopy.gro -gro_frames 0 0 0 0 0 0 0 0 0 0 
#ls -lh $TEST
#vmd -f local_traj.gro -f $RMSD -e loadrmsd.tk


#for i in {0..19}; do
#    for j in {0..9}; do
#        python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i.$j -rmsd_scaletime -rmsd_type angle
#        mv default_$i."$j"ps.png image/2014-09-16-angle
#    done
#done

for i in ${times[@]}; do
    echo $i
done

#for i in ${times[@]}; do
#    python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_scaletime -rmsd_type angle -rmsd_colorrange .7 1.7 -rmsd_grid_file narrow
#    mv narrow.png image/2014-09-18-angle/narrow_"$i"ps.png
#done

#for i in ${times[@]}; do
#    python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_scaletime -rmsd_type angle -rmsd_colorrange 0.0 1.7 -rmsd_grid_file wide
#    mv wide.png image/2014-09-18-angle/wide_"$i"ps.png
#done


#times=(2.5 3.0 3.5 4.0 5.0 7.5 10.0)
#
#DIR=image/2014-09-24-activity
#mkdir $DIR
#for i in ${times[@]}; do
#    FNAME=angle_act
#    python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_type angle -rmsd_activityparm 1.2 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done
#
#times=(0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0 5.0 7.5 10.0)
#
#for i in ${times[@]}; do
#    FNAME=rmsd_act
#    python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_type rmsd -rmsd_activityparm .38 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done
#
#for i in ${times[@]}; do
#    FNAME=diff_act
#    python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_type rmsd -rmsd_scaletime -rmsd_activityparm .12 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
#    mv $FNAME.png $DIR/$FNAME."$i"ps.png
#done


times=(19.9)

for i in ${times[@]}; do
    FNAME=rmsd_act
    #python RMSDAnalyze.py -rmsd_grid png -rmsd_dt $i -rmsd_type rmsd -rmsd_activityparm .38 20 -rmsd_grid_file $FNAME -rmsd_colorrange .1 .9
    python main.py -rmsd_grid display -rmsd_dt $i -rmsd_type angle -rmsd_activityparm 1.2 20 -rmsd_colorrange .1 .9
    #python main.py -op_grid display -op_type q6 -op_grid_file $FNAME 
    #mv $FNAME.png $DIR/$FNAME."$i"ps.png
done

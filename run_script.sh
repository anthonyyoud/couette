#!/bin/sh
#$ -S /bin/sh
#$ -j y -o run$JOB_ID.log
#$ -cwd
#DATA_DIR=$HOME/giga_runs/anom_mag/tmp
DATA_DIR=`pwd`
RUN_DIR=/work/$JOB_ID
#touch $DATA_DIR/RUNNING
mkdir $RUN_DIR
cp couette_mod.out parameters.f90 $RUN_DIR
if [ -e $DATA_DIR/"end_state.dat" ]
then
mv end_state.dat $RUN_DIR
fi
cd $RUN_DIR
mkdir $DATA_DIR/$JOB_ID
ln -s $RUN_DIR/u_growth.dat $DATA_DIR/$JOB_ID/tmp_growth.dat
time ./couette_mod.out
#sleep 10
#if [ ! -e $DATA_DIR/"RUNNING" ]
#then
cp -r $RUN_DIR $DATA_DIR/
mv $DATA_DIR/run$JOB_ID.log $DATA_DIR/$JOB_ID/
rm $DATA_DIR/$JOB_ID/tmp_growth.dat
rm -rf $RUN_DIR
#fi

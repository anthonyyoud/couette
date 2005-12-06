#!/bin/sh
#$ -S /bin/sh
#$ -j y -o run$JOB_ID.log
#$ -cwd

#****************************************************************************
#PARAMETERS
#****************************************************************************
SERIAL=0
NPROCS=4
EXE=couette_mod.out
DATA=parameters.f90
NODECPU=2
NODELIST= 'giga20 giga21' #giga20 giga21' #giga03 giga04 giga06 giga07' #giga10'
#          giga11 giga12 giga13 giga14 giga15 giga16 giga17'
#NODELIST='giga02'
DATA_DIR=`pwd`
RUN_DIR=/work/ay_$JOB_ID

#****************************************************************************
#NO CHANGES NECESSARY BELOW HERE
#****************************************************************************

mkdir $RUN_DIR
mv $EXE $RUN_DIR
cp $DATA $RUN_DIR
if [ -e $DATA_DIR/end_state.dat ]; then
   mv end_state.dat $RUN_DIR
fi
cd $RUN_DIR
mkdir $DATA_DIR/$JOB_ID
if [ $SERIAL == 1 ]; then
   time ./$EXE
else
   echo No. processors = $NPROCS
   echo `hostname` cpu=$NODECPU >> hostfile
   for NODE in $NODELIST
   do
      if [ `hostname` != $NODE.ncl.ac.uk ]; then
         echo $NODE.ncl.ac.uk cpu=$NODECPU >> hostfile
         ssh $NODE "mkdir $RUN_DIR"
         scp $EXE $DATA $NODE:$RUN_DIR
      fi
   done
   lamboot -v hostfile
   lamnodes > nodes.dat
   time mpirun -np $NPROCS $EXE
   lamhalt
fi
rm $EXE
bzip2 -9 `find $RUN_DIR/ -type f ! -iname "end_state.dat" \
                                 ! -iname "parameters.f90" \
                                 ! -iname "xsect*.dat"`

NXSECT=`ls | grep "xsect" | wc -l`
if [ $NXSECT -gt 1 ]; then
  tar cvvf xsect.tar xsect*
  bzip2 -9 xsect.tar
  rm xsect*.dat
elif [ $NXSECT == 1 ]; then
  bzip2 -9 xsect*.dat
fi

cp -r $RUN_DIR/* $DATA_DIR/$JOB_ID/
mv $DATA_DIR/run$JOB_ID.log $DATA_DIR/$JOB_ID/
if [ $SERIAL == 0 ]; then
   for NODE in $NODELIST
   do
      ssh $NODE "cp -r $RUN_DIR $DATA_DIR/;
                 rm -rf $RUN_DIR"
   done
fi
rm -rf $RUN_DIR

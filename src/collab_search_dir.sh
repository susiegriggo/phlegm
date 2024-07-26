#!/bin/bash
#PBS -P ob80
#PBS -lncpus=12
#PBS -lmem=100GB
#PBS -lwalltime=3:00:00
#PBS -ljobfs=1GB
#PBS -l wd
#PBS -l other=hyperthread
#PBS -l storage=gdata/ob80+gdata/if89
#PBS -N collabfold
 
module use /g/data/if89/apps/modulefiles
module load colabfold_batch/1.5.2
 

if [ -z $DIR ]; then
	echo "Please define the DIR environment variable with the -v option" >&2;
	echo "qsub -v DIR=fasta collab_search_dir.sh" >&2;
	exit;
fi

if  [ ! -e $DIR ]; then
	echo "$DIR was not found. Please check" >&2;
	exit;
fi


for FAFILE in $(find $DIR -type f); do
	OUT=`basename $FAFILE .fasta`;
	OUT=$(echo $OUT | sed -e 's/.faa//');
	echo "Output will be in $OUT" >&2;


	# colabfold_search --db-load-mode 1 --threads 12 $FAFILE $COLABFOLDDIR/database $OUT
	qsub -v OUT=$OUT,FAFILE=$FAFILE /home/584/sg5247/pbs/collab_search_only.sh

done

#!/bin/bash -e
#
# *** SUBMITS FML TO CLUSTER***
# *** submits chunks of experiments to cluster ***
# *** Created 27/07/2017 by JDW ***

cd ~/scripts/TMM_model/experiments

#filenames=(20170822_*) # make array with filenames
filenames=($1) # make array with filenames

cd ~/scripts/TMM_model/code
#for i in {$2..$3}
for ((i=$2; i<=$3;i++));
do
	#echo $eachfile
	echo 'submitting:' ${filenames[$i]}
	qsub -j y -o ~/scripts/TMM_model/log -V -S /bin/bash ../code/submit_FML.sh ${filenames[$i]} 3000 20170727_lhs_global_1.0
	sleep 1
done
qstat





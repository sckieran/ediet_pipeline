#!/bin/bash

dir=$1
cutoff=$2
gene=$3
max_jobs=$4
user=$5
minlen=$6
env_name=$7

cat ${dir}/slurm_template.txt ${dir}/scripts/run_seqs.sh > ${dir}/scripts/run_seqs_full.sh

cd ${dir}/${gene}
mkdir -p ./unfiltered_seqfiles
mv fx_col.*.err ./err_and_outs/
mv fx_col.*.out ./err_and_outs/

ls *_clustered.fasta > collapselist
num_seqs=$( wc -l collapselist | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }')
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to make seqfiles for and $tot_per_file sample(s) per job."

split -n l/${max_jobs} --numeric-suffixes=1 collapselist collapselist_
rm collapselist

ls *_filtered_seqs.txt > outslist

num_outs=$( wc -l outslist | awk '{print $1}')
num_jobs=$(ls mksq.*.err | wc -l | awk '{print $1}')

while [[ $num_seqs -ne $num_outs && $num_jobs -le 1000 ]];
do
  for fil in collapselist_*;
  do
    while read p;
    do
      base=$( echo $p | awk -F"_clustered.fasta" '{print $1}')
      if [[ ! -f ${base}_filtered_seqs.txt ]];
      then
        echo "$p" >> temp_$fil
		echo "$p not done"
      fi
    done < ${fil}
    if [[  -s temp_${fil} ]];
    then
    mv temp_${fil} ${fil}
      while true;
	  do
   		echo "outfile for at least one sample in $fil does not yet exist or is empty. Doing $fil."
	 	res=$(sbatch ${dir}/scripts/run_seqs_full.sh $fil ${dir} ${gene} ${cutoff} ${minlen} ${env_name})
   		if squeue -u $user | grep -q "${res##* }"; 
   		then
   		  echo "job ${res##* } for $fil submitted successfully."
       			break
     	elif [[ -f mksq.${res##* }.err ]];
	  	then
	  		echo "job ${res##* } for $fil submitted successfully."
     		break
    	else
	 		echo "job ${res##* } did not submit. Trying again."
	 	fi
  	  done
    else
       echo "all samples for $fil completed already."
    fi
  done
  while true;
  do
        sleep 3s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "mksq" | wc -l | awk '{print $1}')
	echo "waiting for jobs to finish. There are $check jobs remaining."
        if [[ $check -eq 0 ]];then
           echo "done with collapsing ASVs" 
           break
        fi 
  done
  ls *_filtered_seqs.txt > outslist
  num_outs=$( wc -l outslist | awk '{print $1}')
  num_jobs=$(ls mksq.*.err | wc -l | awk '{print $1}')
  echo "there are $num_seqs sequences to make seqfiles for and $num_outs seqfiles. If these numbers don't match, will resubmit jobs as necessary. If these numbers do match, moving on to BLASTing your sequences."
done
mv mksq*.err ./err_and_outs/
mv mksq*.out ./err_and_outs/
mv *_clustered.fasta ./collapsed/
mv *_filtered_seqs.txt ./seqfiles/
rm collapselist_*

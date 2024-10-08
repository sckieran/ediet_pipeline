#!/bin/bash


dir=$1
pattern=$2
r2_pattern=$3
max_jobs=$4
user=$5
gene=$6
env_name=$7

cat ${dir}/slurm_template.txt ${dir}/scripts/run_collapser.sh > ${dir}/scripts/run_collapser_full.sh

cd ${dir}/${gene}

#clean up files#
mkdir -p unpaired paired collapsed seqfiles paired/unassembled  paired/discarded err_and_outs
mv *${pattern} ./unpaired/
mv *${r2_pattern} ./unpaired/
mv *.discarded.fastq ./paired/discarded/
mv *.unassembled*.fastq ./paired/unassembled/
mv pears.*.err ./err_and_outs/
mv pears.*.out ./err_and_outs
rm seqlist_* 
rm outslist

#make list of files to collapse#
ls *_paired.assembled.fastq > pairedlist
num_seqs=$( wc -l pairedlist | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }')
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to cluster and $tot_per_file sample(s) per job."

#cut into slurm jobs for faster processing#
x=1
while [[ $x -le ${max_jobs} ]];
do
  if [[ -s pairedlist ]];
  then
    head -n ${tot_per_file} pairedlist > pairedlist_${x}
    sed -i "1,${tot_per_file}d" pairedlist
    x=$(( $x + 1 ))
  else
    x=$(( $max_jobs + 1 ))
  fi
done
rm pairedlist
while [[ $num_seqs -ne $num_outs ]];
do
  for fil in pairedlist_*;
  do
    while read p;
    do
      base=$( echo $p | awk -F"_paired.assembled.fastq" '{print $1}')
      if [[ ! -f ${base}_clustered.fasta ]];
      then
        echo "$p" >> temp_$fil
      fi
    done < ${fil}
    if [[  -s temp_${fil} ]];
    mv temp_${fil} ${fil}
    then
      while true;
     	do
     		echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     		res=$(sbatch ${dir}/scripts/run_collapser_full.sh $fil ${dir}/${gene} ${env_name})
   		if squeue -u $user | grep -q "${res##* }"; 
   		then
   			echo "job ${res##* } for $fil submitted successfully."
    			break
     		elif [[ -f fx_col.${res##* }.err ]];
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
        check=$(echo "$chck" | grep "fx_col" | wc -l | awk '{print $1}')
        if [[ $check -eq 0 ]];then
           echo "done with collapsing ASVs" 
           break
        fi 
  done
  ls *_clustered.fasta > outslist
  num_outs=$( wc -l outslist | awk '{print $1}')
done

echo "there are $num_seqs input samples and $num_outs clustered output samples. Moving on."

rm outslist
rm pairedlist_*
mv *_paired.assembled.fastq ./paired/
rm fx_col.*.err
rm fx_col.*.out

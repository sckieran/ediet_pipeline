#!/bin/bash

pattern=$1
r2_pattern=$2
dir=$3
max_jobs=$4
user=$5
gene=$6
env_name=$7

cd ${dir}/${gene}

cat ${dir}/slurm_template.txt ${dir}/scripts/pear.sh > ${dir}/scripts/pear_full.sh

 
ls *${pattern} > seqlist
num_seqs=$( cat seqlist | wc -l | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }' )
echo "max_jobs is set to $max_jobs"
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to pear and $tot_per_file sample(s) per job."

x=1
while [[ $x -le ${max_jobs} ]];
do
 # echo "x is $x and max_jobs is $max_jobs"
  if [[ -s seqlist ]];
  then
   # echo "seqlist not empty, making seqlist_${x}"
    head -n ${tot_per_file} seqlist > seqlist_${x}
    sed -i "1,${tot_per_file}d" seqlist
    rem=$( wc -l seqlist | awk '{print $1}')
   # echo "there are $rem samples remaining in the seqlist"
    x=$(( $x + 1 ))
  else
    echo "seqlist empty, moving on"
    x=$(( ${max_jobs} + 1 ))
  fi
done
rm seqlist
num_outs=1
while [[ $num_outs -ne $num_seqs ]];
do
  for fil in seqlist_*;
  do  
    while read p;
    do
      base=$( echo $p | awk -F"${pattern}" '{print $1}')
      if [[ ! -f ${base}_paired.assembled.fastq ]];
      then 
        echo "$p" >> temp_$fil
      fi
    done < ${fil}
    if [[  -s temp_${fil} ]]
    then
     mv temp_${fil} ${fil}
      while true;
     			do
     				echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     				res=$(sbatch ${dir}/scripts/pear_full.sh $fil $pattern $r2_pattern ${dir}/${gene} ${env_name})
   			  	if squeue -u $user | grep -q "${res##* }"; 
   	  			then
     					echo "job ${res##* } for $fil submitted successfully."
       					break
       				elif [[ -f pears.${res##* }.err ]];
	  			then
	  			  	echo "job ${res##* } for $fil submitted successfully."
      					break
      				else
	  				echo "job ${res##* } did not submit. Trying again."
				fi
   			done
    else
      echo "all files for $fil completed."
    fi
  done
  while true;
  do
       sleep 3s
       ck="squeue -u ${user}"
       chck=$($ck)
       check=$(echo "$chck" | grep "pear" | wc -l | awk '{print $1}')
       if [[ $check -eq 0 ]];then
          echo "done with pears" 
          break
       fi 
   done
  ls *_paired.assembled.fastq > outslist
  num_outs=$( wc -l outslist | awk '{print $1}')
done

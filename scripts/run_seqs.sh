#!/bin/bash

#SBATCH -J mksq
#SBATCH -e mksq.%j.out
#SBATCH -o mksq.%j.err

inp=$1
dir=$2
gene=$3
rlib=$4
cutoff=$5
minlen=$6
env_name=$7

cd ${dir}/${gene}


while read fil;
do
  base=$(echo $fil | awk -F"_clustered.fasta" '{print $1}')
  grep ">" $fil | awk -F"-" '{print $2}' > temp_reads_${base}
  grep -v ">" $fil  > temp_seqs_${base}
  paste  temp_seqs_${base} temp_reads_${base}  > ${base}_seqs.txt
  cat ${base}_seqs.txt | awk -v m=$minlen '{ if (length($1) > m) print }' > temp_${base}_seqs.txt
  mv temp_${base}_seqs.txt ${base}_seqs.txt
  eval "$(conda shell.bash hook)"
	conda activate $env_name
  Rscript ${dir}/scripts/filter_rra.R ${dir}/${gene} ${base}_seqs.txt $base $cutoff
  mv ${base}_seqs.txt ./unfiltered_seqfiles/${base}_seqs.txt
  rm temp_reads_${base} temp_seqs_${base}
done < ${inp}


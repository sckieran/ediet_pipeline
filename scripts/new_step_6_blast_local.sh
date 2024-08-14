#!/bin/bash

while getopts ":n:g:d:m:r:b:c:t:e:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) gene="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    m) minlen="$OPTARG"
    ;;
    r) db_dirr="$OPTARG"
    ;;
    b) localdat="$OPTARG"
    ;;
    c) cutoff="$OPTARG"
    ;;
    t) return_low="$OPTARG"
    ;;
    e) env_name="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    esac
done

##load conda
conda_source=$(conda info | grep -i 'base environment' | awk -F":" '{print $2}'| awk -F" " '{print $1}')
source ${conda_source}/etc/profile.d/conda.sh
conda activate ${env_name}
#source activate ${env_name}
pref=$CONDA_PREFIX

cd ${dirr}

if [[ -z ${localdat} ]]
then
	localdat=${prefix}_${gene}_reference
fi

mkdir -p ${gene}_out

cd ${gene}_out
cp ${dirr}/${gene}/seqfiles/* .

echo "copied files. Beginning blast and taxtable."


cd ${dirr}/${gene}_out
cat *_seqs.txt | cut -f1 | sort | uniq | awk -v m=$minlen '{ if (length($0) > m) print }' > temp_seqs
sed -i '/^$/d' temp_seqs

echo "making query fasta from sample sequences"
	##make query fasta from seqlist#
x=1
n=$(wc -l temp_seqs | awk '{print $1}')
touch ${prefix}_${gene}_headers
while [[ $x -le $n ]]
do
	if [[ $x -le 9 ]]
 	then 
		echo ">seq_00000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 99 ]] && [[ $x -ge 10 ]]
	then
		echo ">seq_0000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 999 ]] && [[ $x -ge 100 ]]
	then 
		echo ">seq_000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 9999 ]] && [[ $x -ge 1000 ]]
	then
		echo ">seq_00${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 99999 ]] && [[ $x -ge 10000 ]]
	then
		echo ">seq_0${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -ge 100000 ]]
	then
		echo ">seq_${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	fi
done
paste -d '\n' ${prefix}_${gene}_headers temp_seqs > ${prefix}_${gene}_combined_ASVs.fasta
rm temp_seqs ${prefix}_${gene}_headers

echo "done with query sequences, looking for NCBI taxonomy file"
cd ${dirr}
if ls ncbi*.csv* 1> /dev/null 2>&1; then
	echo "ncbi tax file found, beginning tax assessment"
	cp ncbi*.csv* ${dirr}/${gene}_out/
	cd ${dirr}/${gene}_out/
	gunzip ncbi*.csv.gz
else
	echo "no ncbi tax file found, attempting to install and run ncbitax2lin."
	pip install -U ncbitax2lin
	wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump
	ncbitax2lin --nodes-file taxdump/nodes.dmp --names-file taxdump/names.dmp
	cp ncbi*.csv* ${dirr}/${gene}_out/
	cd ${dirr}/${gene}_out/
	gunzip ncbi*.csv.gz
fi

cd ${dirr}/${gene}_out
ncbi=$( ls ncbi*.csv | head -n1 | awk '{print $1}')
echo "now doing local blast search, this may take many hours. You have $n sequences to align. You can check your progress in the ${gene}_out folder with the command: tail ${prefix}_${gene}_raw_blast_out"
echo "localdat is $localdat, prefix is $prefix"

if [[ ${return_low} == "TRUE" ]]
then
	echo "you set return_low to TRUE, so BLAST will return the top bitscore matches regardless of percent identity."
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 10 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out
else
 	echo "you set return_low to FALSE, or did not enter a valid TRUE/FALSE value, so BLAST will only return hits above ${cutoff} percent identity, regardless of score."
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 10 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out -perc_identity ${cutoff}
  fi
blastout=${prefix}_${gene}_raw_blast_out
sed -i '/^#/d' $blastout
echo "blast done, making outfiles" 
 
for f in *_seqs.txt; do awk "{print \"$f\t\" \$0}" "$f"; done > cat_file_list.txt

sed -i 's/_seqs.txt//g' cat_file_list.txt
sed -i 's/_filtered//g' cat_file_list.txt

paste - - < ${prefix}_${gene}_combined_ASVs.fasta > asvs.txt

Rscript ${dirr}/scripts/get_best_hits.R ${prefix} ${gene} ${pref}/lib/R/library/ ${dirr}/${gene}_out/ ${ncbi}

mkdir -p ${dirr}/${gene}_out/sample_seqfiles ${dirr}/results_tables
mv *_seqs.txt ${dirr}/${gene}_out/sample_seqfiles
cd ${dirr}
cp ${dirr}/${gene}_out/*_taxatable.txt ${dirr}/results_tables/

#!/bin/bash

while getopts ":n:g:d:m:r:b:c:t:j:u:" opt; do
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
    j) max_jobs="$OPTARG"
    ;;
    u) user="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    esac
done


module load ncbi-blast


cd ${dirr}

if [[ -z ${localdat} ]]
then
	localdat=${prefix}_${gene}_reference
fi

mkdir -p ${gene}_out 

cd ${gene}_out
cp ${dirr}/${gene}/seqfiles/* .
 	
echo "copied files. Beginning blast and taxtable."
cat *_seqs.txt | cut -f1 | sort | uniq > temp_seqs
sed -i '/^$/d' temp_seqs

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
	elif [[ $x -le 999999 ]]
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
cd ${dirr}

##get or make taxfile#
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
echo "now doing local blast search, this may some time. You have $n sequences to align. You can check your progress in the ${gene}_out folder with the command: tail ${prefix}_${gene}_raw_blast_out"
echo "localdat is $localdat, prefix is $prefix"

if [[ ${return_low} == "TRUE" ]]
then
	echo "you set return_low to TRUE, so BLAST will return the top bitscore matches regardless of percent identity."
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 3 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out
else
 	echo "you set return_low to FALSE, or did not enter a valid TRUE/FALSE value, so BLAST will only return hits above ${cutoff} percent identity, regardless of score."		
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 3 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out -perc_identity ${cutoff}
  fi
blastout=${prefix}_${gene}_raw_blast_out
sed -i '/^#/d' $blastout
echo "blast done, making outfiles"

#make a list of no-hits#
grep ">" ${prefix}_${gene}_combined_ASVs.fasta | awk -F">" '{print $2}' | sort > out1
cut -f1 ${blastout}| sort | uniq > out2
comm -23 out1 out2 > list_of_no_hits
totalseqs=$( wc -l out1 | awk '{print $1}')
totalhits=$( wc -l out2 | awk '{print $1}')
nohits=$( wc -l list_of_no_hits | awk '{print $1}')
echo "there were ${totalhits} raw BLAST hits out of ${totalseqs} unique sequences, and ${nohits} raw BLAST no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."

if [[ ${totalhits} -eq 0 ]]
then
	echo "there were no local blast hits. This is a problem. Exiting. Check your query fasta, ${prefix}_${gene}_combined_ASVs.fasta and your seqs.txt files for errors."
     	exit 1;
fi
echo "done with local blast, now doing remote blast. There are ${n} sequences to align. This may take many hours. This option is not recommended if you have >50,000 sequences to align. Do NOT set taxa_rra to 0 and then choose this option."
if [[ ${return_low} == "TRUE" ]]
then
	echo "you set return_low to TRUE, so BLAST will return the top bitscore matches regardless of percent identity."
  	blastn -db nt -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore staxids" -culling_limit 1 -out ${prefix}_${gene}_remote_raw_blast_out -remote
else
 	echo "you set return_low to FALSE, or did not enter a valid TRUE/FALSE value, so BLAST will only return hits above ${cutoff} percent identity, regardless of score."		
  	blastn -db nt -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore staxids" -culling_limit 1 -out ${prefix}_${gene}_remote_raw_blast_out -perc_identity ${cutoff} -remote
  fi
#make a list of your unique sequences

remote_blastout=${prefix}_${gene}_remote_raw_blast_out
sed -i '/^#/d' $blastout
echo "blast done, making outfiles"

#make a list of no-hits#
cut -f1 ${remote_blastout}| sort | uniq > out2
comm -23 out1 out2 > remote_list_of_no_hits
totalseqs=$( wc -l out1 | awk '{print $1}')
totalhits=$( wc -l out2 | awk '{print $1}')
nohits=$( wc -l remote_list_of_no_hits | awk '{print $1}')
echo "there were ${totalhits} raw BLAST hits out of ${totalseqs} unique sequences, and ${nohits} remote BLAST no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."
cp out1 temp_seqlist
rm out1 out2

##add your species and taxid to your hits##
echo "adding species and taxid to your remote blast results"
cut -f5 $remote_blastout | awk -F" " '{print $1,$2}' > temp_spec
cut -f7 $remote_blastout > temp_taxids
cut -f6 $remote_blastout > temp_scores
cut -f1-4 $remote_blastout | paste - temp_spec temp_taxids temp_scores > ${remote_blastout}_with_tax
rm temp_spec temp_taxids temp_scores

echo "adding species and taxid to your local blast results"
cut -f5 $blastout | awk -F" " '{print $1,$2}' > temp_spec
cut -f5 $blastout | awk -F"taxid=" '{print $2}' | awk -F" " '{print $1}' > temp_taxids
cut -f6 $blastout > temp_scores
cut -f1-4 $blastout | paste - temp_spec temp_taxids temp_scores > ${blastout}_with_tax
rm temp_spec temp_taxids temp_scores

##modify your ncbi tax file to contain only taxa within your reference database##
echo "now modifying your taxonomy file to limit your search space. This saves time."
cut -f6 ${blastout}_with_tax | sort | uniq > all_local_taxids
cut -f6 ${remote_blastout}_with_tax | sort | uniq > all_remote_taxids
cat all_local_taxids all_remote_taxids | sort | uniq > all_taxids
sed -i 's/$/,/g' all_taxids
sed -i 's/^/\^/g' all_taxids
grep -f all_taxids $ncbi | cut -f1,3,4,5,6,7,8 -d"," > ${ncbi}_r
sed -i 's/,/\t/g' ${ncbi}_r
rm all_remote_taxids all_local_taxids all_taxids
tot=$( wc -l temp_seqlist | awk '{print $1}')

 echo "making your taxonomic assignment files and beginning to assign jobs."
tot_per_file=$( awk -v a1=$tot -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }' )
if [[ ${tot_per_file} -eq 0 ]]
	then
 		 tot_per_file=1
	fi
echo "there were $tot samples to assign taxonomy for and $tot_per_file sample(s) per job."
x=1
while [[ $x -le ${max_jobs} ]];
do
	if [[ -s temp_seqlist ]];
  	then
   		head -n ${tot_per_file} temp_seqlist > seqlist_${x}
    		sed -i "1,${tot_per_file}d" temp_seqlist
   		x=$(( $x + 1 ))
 	else
  		x=$(( $max_jobs + 1 ))
  	fi
done
rm temp_seqlist
tot2=$(( $tot * 2 ))
while [[ $num_outs -ne $tot2 ]];
do
	for fil in seqlist_*;
 	do
  		x=$( echo $fil | awk -F"_" '{print $2}')
    		while true;
     		do
     			if [[ ! -s ${prefix}_${gene}_best_blast_hits.out_${x} ]];
	 		then
	 			echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     				res=$(sbatch ${dirr}/scripts/run_tax.sh $x $prefix $gene $tot_per_file $blastout $ncbi $dirr $remote_blastout)
   				if [[ ! -s remote_${prefix}_${gene}_best_blast_hits.out_${x} ]];
       				then
       					sbatch ${dirr}/scripts/run_tax_remote.sh $x $prefix $gene $tot_per_file $remote_blastout $ncbi $dirr
       				fi
	   			if squeue -u $user | grep -q "${res##* }"; 
   				then
   					echo "job ${res##* } for $fil submitted successfully."
       					break
       				elif [[ -f tax.${res##* }.err ]];
	  			then
	  				echo "job ${res##* } for $fil submitted successfully."
      					break
      				else
	  				echo "job ${res##* } did not submit. Trying again."
				fi
     			elif [[ ! -s remote_${prefix}_${gene}_best_blast_hits.out_${x} ]];
			then
   				sbatch ${dirr}/scripts/run_tax_remote.sh $x $prefix $gene $tot_per_file $blastout $ncbi $dirr
       				break
			else
	 			echo "all samples from $fil already done."
     				break
      			fi
    		done
       	done	
	while true;
     	do
       		sleep 3s
	 	ck="squeue -u ${user}"
		chck=$($ck)
  		check=$(echo "$chck" | grep "tax" | wc -l | awk '{print $1}')
		echo "waiting for jobs to run. There are $check jobs left"
       		if [[ $check -eq 0 ]];
	 	then
        		echo "no jobs left, checking that jobs ran successfully." 
         		break
       		fi 
	done
  	cat *_best_blast_hits.out_* > outslist
   	num_outs=$( wc -l outslist | awk '{print $1}')
    	if [[ $num_outs -gt $tot2 ]]
     	then
      		num_outs=$tot2
	fi
    	echo "there are $tot2 sequences to assign and $num_outs sequences successfully assigned. If these numbers match, moving on to taxonomy. If not, checking each sample and re-submitting jobs as needed."
done
    
#cat your files and make a header for the best hits table.

for fil in ${prefix}_${gene}_best_blast_hits.out_*;
do
 	cut -f3-13 remote_${fil} > temp_remote
  	paste $fil temp_remote > temp_fil
  	mv temp_fil $fil
   	rm temp_fil temp_remote
done


cat ${prefix}_${gene}_best_blast_hits.out_* | sort -k1 > ${prefix}_${gene}_best_blast_hits.out
echo "done with choosing best blast hits, now creating and formatting outfiles."
	
echo "sequence	seqnum	local_identity	local_species	local_taxid	local_phylum	local_class	local_order	local_family	local_genus	local_bitscore	local_num_spec_in_best_hit	local_all_spec_in_best_hit	remote_identity	remote_species	remote_taxid	remote_phylum	remote_class	remote_order	remote_family	remote_genus	remote_bitscore	remote_num_spec_in_best_hit	remote_all_spec_in_best_hit" > ${prefix}_${gene}_best_blast_hits.header
cat ${prefix}_${gene}_best_blast_hits.header ${prefix}_${gene}_best_blast_hits.out > ${prefix}_${gene}_best_blast_hits.txt

#clean up outfiles
rm ${prefix}_${gene}_best_blast_hits.out*
rm seqlist_*
rm tax.*.err
rm tax.*.out

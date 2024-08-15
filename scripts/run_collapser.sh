#SBATCH -J fx_col
#SBATCH -e fx_col.%j.err
#SBATCH -o fx_col.%j.out

infil=$1
dir=$2
env_name=$3

##load conda
conda_source=$(conda info | grep -i 'base environment' | awk -F":" '{print $2}'| awk -F" " '{print $1}')
source ${conda_source}/etc/profile.d/conda.sh
conda activate ${env_name}
#source activate ${env_name}



cd $dir

while read p;
do
  base=$( echo $p | awk -F"_paired.assembled.fastq" '{print $1}')
  fastx_collapser -v -i $p -o ${base}_clustered #collapse ASVs
  doub_num=$( grep -m1 -n ">*-2$" ${base}_clustered | awk -F":" '{print $1}') #remove singletons and doubletons
  head_num=$(( $doub_num - 1 ))
  head -n $head_num ${base}_clustered > ${base}_clustered.fasta #rename files
  rm ${base}_clustered
done < $infil

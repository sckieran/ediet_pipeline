#SBATCH -J pear_j
#SBATCH -e pears.%j.err
#SBATCH -o pears.%j.out


infil=$1
pattern=$2
r2_pattern=$3
dir=$4
env_name=$5

cd ${dir}
##load conda
conda_source=$(conda info | grep -i 'base environment' | awk -F":" '{print $2}'| awk -F" " '{print $1}')
source ${conda_source}/etc/profile.d/conda.sh
conda activate ${env_name}
#source activate ${env_name}


while read fil;
do
  base=$(echo $fil | awk -F"$pattern" '{print $1}')
  pear -f ${fil} -r ${base}${r2_pattern} -o ${base}_paired -j 10
done < $infil

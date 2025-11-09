#Tisha Melia
#This file for assembling the chloroplast genome
#!/bin/bash

#SBATCH --job-name=trim_galore
#SBATCH --ntasks=5
##SBATCH --partition=short
#SBATCH --cpus-per-task=4


############################################################################
# SETTING VARS
############################################################################

# Now let's keep track of some information just in case anything goes wrong
start=$(date +%s%N)
echo "=========================================================="
echo "Starting on : ${start}"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "=========================================================="

PARALLEL=4
   
HOME="/mgpfs/home/tishalia/"

module load bioinformatics/anaconda3
module load bioinformatics/blast

eval "$(conda shell.bash hook)"
conda activate getorganelle

outdir=${HOME}/data/NGS-509-WGS-Illumina
read1=${outdir}/N509-1_S8_R1_001_val_1.fq.gz
read2=${outdir}/N509-1_S8_R2_001_val_2.fq.gz


echo "get_organelle_from_reads.py -1 $read1 -2 $read2 -o ${HOME}/assembly/default --max-reads 9E7 -t 20 -R 15 -k 21,45,65,85,105 -F embplant_pt --overwrite"
get_organelle_from_reads.py -1 $read1 -2 $read2 -o ${HOME}/assembly/default --max-reads 9E7 -t 20 -R 15 -k 21,45,65,85,105 -F embplant_pt --overwrite


echo "done running"
echo "=========================================================="

# print out some diagnostic stuff
end=$(date +%s%N)
duration=$(((end - start) / 1000000000))

echo "Stop time is ${end}"
echo "Duration: ${duration}"



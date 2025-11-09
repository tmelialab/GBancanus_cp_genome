#Tisha Melia
#Script to perform qc with fastqc

#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
##SBATCH --partition=short
#SBATCH --cpus-per-task=4
#
# Tisha Melia
# May 24, 2024: 
# Do fastqc result

############################################################################
# SETTING VARS
############################################################################

# Now let's keep track of some information just in case anything goes wrong
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "=========================================================="


#Get command line arguments
INPUT_DIR=$1
OUTPUT_DIR=$2

echo $INPUT_DIR
echo $OUTPIT_DIR
cd $OUTPUT_DIR
FASTQC="/mgpfs/home/tishalia/tool/FastQC/fastqc"


# print out some diagnostic stuff
echo "Running on host $(hostname)"
echo ""
echo "Directory is $(pwd)"
echo ""
echo "Start time is $(date)"
echo "Start time is $(date)" > $LOG
echo ""


# run my commands
echo ""
echo "running commands..."
echo ""


for f in $INPUT_DIR/*.gz
do
  echo "Processing $f file..."
  echo "$FASTQC $f"
  $FASTQC $f
done



echo "done running"
echo "=========================================================="
#Just list everything in scratch dir (Don't know what output files I need)
#Also want the file size of everything
echo "Size of files in scratch"
du -ch *

# print out some diagnostic stuff
echo "Stop time is $(date)"
echo "Stop time is $(date)" >> $LOG
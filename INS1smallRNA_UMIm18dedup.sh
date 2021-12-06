#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J UMIdedup
#SBATCH --output=UMIm18dedup-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=128g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-5

module load umi-tools/1.0.0
 
names=($(cat smallRNAexpList_2019.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:53}
INPUT_BAM=/bgfs/rnicholls/INS1_PWS_smallRNA/Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}_sortedbycoord.bam

echo $READ_BASE
echo $INPUT_BAM

umi_tools dedup \
-I $INPUT_BAM \
-L ./Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}m18.log \
--method percentile \
-S ./Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}_UMIm18dedup_.bam

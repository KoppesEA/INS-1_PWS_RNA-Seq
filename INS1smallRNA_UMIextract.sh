#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J UMIextract
#SBATCH --output=UMIextract-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=48g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-5

module load umi-tools/0.5.5
 
names=($(cat smallRNAexpList_2019))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:53}
INPUT_FASTQ_R1=./smallRNA_raw_fastq/${READ_BASE}_R1_001.fastq.gz

echo $READ_BASE
echo $INPUT_FASTQ_R1

umi_tools extract \
--stdin=$INPUT_FASTQ_R1 \
--log=./UMIextract/UMIsmallRNA_${READ_BASE}.log \
--stdout=./UMIextract/${READ_BASE}_R1_UMIextract.fa.gz \
--extract-method=regex \
--bc-pattern=".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)"


#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J UMIcutAdapt
#SBATCH --output=UMIcutAdaptm18-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=48g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-5

module load cutadapt/1.18
 
names=($(cat smallRNAexpList_2019.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:53}
INPUT_FASTQ_R1=./UMIextract/${READ_BASE}_R1_UMIextract.fq.gz

echo $READ_BASE
echo $INPUT_FASTQ_R1

cutadapt \
-j 4 \
-a AACTGTAGGCACCATCAAT \
-o ./UMIcutAdaptm18/${READ_BASE}_R1_UMIcutAdapt_m18.fq.gz \
-q 15,10 \
-m 18 \
$INPUT_FASTQ_R1


#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J RNA_QCTRIM
#SBATCH --output=RNA_QCTRIM-%A_%a.out
##array should start from zero
##SBATCH --array=0-5

module load FastQC/0.11.5
module load TrimGalore/0.4.5

names=($(cat expList_2017.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

INPUT_FASTQ1=${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz
INPUT_FASTQ2=${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz
OUTPUT=${names[${SLURM_ARRAY_TASK_ID}]:60}


echo $INPUT_FASTQ1
echo $INPUT_FASTQ2
echo $OUTPUT

#1. Fastqc passed 
fastqc -o $OUT_DIR/$INPUT_FASTQ1 $INPUT_FASTQ1.fq
fastqc -o $OUT_DIR/$INPUT_FASTQ2 $INPUT_FASTQ2.fq


#2. Trim_galore: cutadapt
mkdir ./$OUTPUT
trim_galore --paired --retain_unpaired --output ./$OUTPUT \
$INPUT_FASTQ1 $INPUT_FASTQ2


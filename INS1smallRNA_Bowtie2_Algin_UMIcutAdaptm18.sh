#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J Bowtie2smallRNA
#SBATCH --output=Bowtie2smallRNAUMIcutadapm18-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-5

module load gcc/8.2.0
module load bowtie2/2.3.4.2
 
names=($(cat smallRNAexpList_2019.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:53}
INPUT_FASTQ_CUTADAPT=/bgfs/rnicholls/INS1_PWS_smallRNA/UMIcutAdaptm18/${READ_BASE}_R1_UMIcutAdapt_m18.fq.gz
BOWTIE_INDEX=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor_v99/Rnor6bt2index/Rnor6bt2index

echo $READ_BASE
echo $INPUT_FASTQ_CUTADAPT
echo $BOWTIE_INDEX

bowtie2 \
-p 4 \
-x $BOWTIE_INDEX \
-U $INPUT_FASTQ_CUTADAPT \
-S ./Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}.sam

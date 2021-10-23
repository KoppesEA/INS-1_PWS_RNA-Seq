#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J STAR
#SBATCH --output=STARcustom_new-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=48g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-5

module load star/2.7.0e
 
names=($(cat expList_2017))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:60}
INPUT_FASTQ1_trim=./$READ_BASE/${READ_BASE}_R1_001_val_1.fq.gz
INPUT_FASTQ2_trim=./$READ_BASE/${READ_BASE}_R2_001_val_2.fq.gz
ANNOTDIR=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor6.0_custom/Rnor_Custom_100519/Rnor6.0_custom_STAR


echo $READ_BASE
echo $INPUT_FASTQ1_trim
echo $INPUT_FASTQ2_trim

# Aligning by STAR
rm -rf ./STARcustom_v98/${READ_BASE}tmp/*
rmdir ./STARcustom_v98/${READ_BASE}tmp
mkdir ./STARcustom_v98/${READ_BASE}
STAR \
--outTmpDir ./STARcustom_v98/${READ_BASE}tmp \
--runThreadN 4 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outSAMprimaryFlag AllBestScore \
--chimSegmentMin 15 \
--chimOutType Junctions \
--twopassMode Basic \
--readFilesIn $INPUT_FASTQ1_trim $INPUT_FASTQ2_trim \
--outFileNamePrefix ./STARcustom_v98/${READ_BASE}/${READ_BASE} \
--quantMode GeneCounts \
--outStd Log \
--genomeDir $ANNOTDIR \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate
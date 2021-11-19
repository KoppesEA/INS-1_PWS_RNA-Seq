#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J RSEMcount
#SBATCH --output=RSEMcount_Rnorcustv98-%A_%a.out
#SBATCH --cpus-per-task=4 # 
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-5

names=($(cat expList_2017))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

WKDIR=/bgfs/rnicholls/hpark_share/EK_100618
READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:60}
INPUT_FASTQ1_trim=$WKDIR/$READ_BASE/${READ_BASE}_R1_001_val_1.fq.gz
INPUT_FASTQ2_trim=$WKDIR/$READ_BASE/${READ_BASE}_R2_001_val_2.fq.gz
ANNOT=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor6.0_custom/Rnor_Custom_100519/Rnor_STAR_RSEM


echo $WKDIR
echo $READ_BASE
echo $INPUT_FASTQ1_trim
echo $INPUT_FASTQ2_trim
echo $ANNOT

module load gcc/8.2.0
module load rsem/1.3.1

rsem-calculate-expression \
	-p 4 \
	--paired-end \
	--strandedness reverse \
	--star \
	--star-path /ihome/crc/install/star/STAR-2.7.3a/bin/Linux_x86_64_static/ \
	--star-gzipped-read-file \
	--sort-bam-by-coordinate \
	$INPUT_FASTQ1_trim \
	$INPUT_FASTQ2_trim \
	$ANNOT \
	$READ_BASE
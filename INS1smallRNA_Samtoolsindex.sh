#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J SamRIP
#SBATCH --output=smallRNA_samtools
#SBATCH --cpus-per-task=4

module load gcc/8.2.0
module load samtools/1.9

WKDIR=/bgfs/rnicholls/INS1_PWS_smallRNA/Bowtie2cutUMIAdaptalignedm18_v99

for bamfile in ${WKDIR}/*.bam
do
echo $bamfile

samtools index $bamfile

done
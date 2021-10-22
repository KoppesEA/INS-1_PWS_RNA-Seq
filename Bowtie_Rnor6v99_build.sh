#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J RnorBowtieBuild 
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem=256g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=RnorBowtieBuild_indexing.out

GENOME=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor_v99/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
OUTPUT=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor_v99/Rnor6bt2index
 
module load gcc/8.2.0
module load bowtie2/2.3.4.2 

bowtie2-build \
--threads 16 \
$GENOME \
$OUTPUT

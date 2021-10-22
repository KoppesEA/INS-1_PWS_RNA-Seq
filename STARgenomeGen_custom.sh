#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J STARcustom
#SBATCH --output=genSTARGenome-%A_%a.out
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem=256g # Memory pool for all cores (see also --mem-per-cpu)

module load star/2.7.0e

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ./Rnor6.0_custom_STAR \
--genomeFastaFiles ./Rattus_norvegicus.Rnor_6.0.dna.custom.fa \
--sjdbGTFfile ./Rattus_norvegicus.Rnor_6.0.98.custom.gtf \
--sjdbOverhang 74

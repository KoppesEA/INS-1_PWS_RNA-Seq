#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J RSEMprep
#SBATCH --output=RSEMprep_Rnorcustv98-%A_%a.out
#SBATCH --cpus-per-task=4 # 
#SBATCH --mem=32g # Memory pool for all cores (see also --mem-per-cpu)


WKDIR=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor6.0_custom/Rnor_Custom_100519
ANNOT=$WKDIR/Rattus_norvegicus.Rnor_6.0.98.custom.gtf
FASTA=$WKDIR/Rattus_norvegicus.Rnor_6.0.dna.custom.fa
STAR=/ihome/crc/install/star/STAR-2.7.3a/bin/Linux_x86_64_static/ #use module load and which STAR to find directory

echo $WKDIR
echo $ANNOT
echo $FASTA

module load gcc/8.2.0
module load rsem/1.3.1

rsem-prepare-reference \
	-p 4 \
	--gtf $ANNOT \
	--star \
	--star-path $STAR \
	$FASTA \
	$WKDIR/Rnor_STAR_RSEM
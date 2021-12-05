#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J htseq-smallRNA1
#SBATCH --output=htseq-smallRNA-cutUMIm18-UnionAllExon-%A_%a.out
#SBATCH --cpus-per-task=1 # HTseq does not support multithread
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-5

names=($(cat smallRNAexpList_2019))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:53}
WKDIR=/bgfs/rnicholls/INS1_PWS_smallRNA/Bowtie2cutUMIAdaptalignedm18_v99
INPUT_BAM=$WKDIR/${READ_BASE}_UMIm18dedup_.bam
ANNOT=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor_v99/Rattus_norvegicus.Rnor_6.0.99.gtf

echo $READ_BASE
echo $INPUT_FASTQ_CUTADAPT
echo $ANNOT

module load htseq/0.11.2

htseq-count \
--format bam \
--order pos \
--mode union \
--nonunique all \
--minaqual 1 \
--stranded yes \
--type exon \
--idattr gene_id \
--additional-attr gene_name \
$INPUT_BAM $ANNOT > $WKDIR/htseq_UnionAll/${READ_BASE}.tsv
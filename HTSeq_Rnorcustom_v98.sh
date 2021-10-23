#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J htseq-Rnorcustv98
#SBATCH --output=htseq_Rnorcustv98-%A_%a.out
#SBATCH --cpus-per-task=1 # HTseq does not support multithread
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-5

names=($(cat expList_2017))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]:60}
INPUT_BAM=/bgfs/rnicholls/hpark_share/EK_100618/STARcustom_v98/$READ_BASE/${READ_BASE}Aligned.sortedByCoord.out.bam
ANNOT=/bgfs/rnicholls/REFGenomes/Rnor_6.0/Rnor6.0_custom/Rnor_Custom_100519/Rattus_norvegicus.Rnor_6.0.98.custom.gtf

echo $READ_BASE
echo $INPUT_BAM
echo $ANNOT

module load htseq/0.11.2
#htseq with union overlap and non-unique none options
htseq-count \
--format bam \
--order pos \
--mode union \
--nonunique all \
--minaqual 1 \
--stranded reverse \
--type exon \
--idattr gene_id \
--additional-attr gene_name \
$INPUT_BAM $ANNOT > ./htseq_custom_unionnone/${READ_BASE}.tsv

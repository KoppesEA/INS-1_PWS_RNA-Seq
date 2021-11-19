#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J RSEMmatrix
#SBATCH --output=RSEMmatrix_Rnorcustv98-%A_%a.out
#SBATCH --cpus-per-task=1

module load gcc/8.2.0
module load rsem/1.3.1

rsem-generate-data-matrix \
Control_1_S1.genes.results \
Control_2_S2.genes.results \
Control_3_S3.genes.results \
PWS_1_S4.genes.results \
PWS_2_S5.genes.results \
PWS_3_S6.genes.results \
> INS1_RSEM_genes.matrix

rsem-generate-data-matrix \
Control_1_S1.isoforms.results \
Control_2_S2.isoforms.results \
Control_3_S3.isoforms.results \
PWS_1_S4.isoforms.results \
PWS_2_S5.isoforms.results \
PWS_3_S6.isoforms.results \
> INS1_RSEM_isoforms.matrix


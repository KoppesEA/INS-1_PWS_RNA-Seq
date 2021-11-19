# INS-1_PWS_PNAS
Markdown on the useage of uploaded bash scripts for RNA-Seq analysis associated with Koppes et al. 2021 for analysis of total- and small- RNA-Seq of control and PWS INS-1 cells.  

Fastq files associated with this study can be found via GEO  
GSExxxx  (Total RNA)  
and  
GSExxx  (small RNA)


*Directories and version numbers may need to be adjusted within each script by user to run as needed*

## Part I: Custom Genome Generation
Download Rnor6 top level genome sequence and gene set annotation from Ensembl (https://useast.ensembl.org/info/data/ftp/index.html)  

    wget http://ftp.ensembl.org/pub/release-98/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
    wget http://ftp.ensembl.org/pub/release-98/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.98.gtf.gz

Note that as of 10/21/2021 current genome and annotation is now Rnor_6.0.104  

Unzip the genome fasta file and then append the custom transgne sequences from `mINS2tg.fa` and `hINStg.fa`

     gunzip -c Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz > Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
     cat Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa hINStg.fa mINS2tg.fa > Rattus_norvegicus.Rnor_6.0.dna.custom.fa
     gunzip -c Rattus_norvegicus.Rnor_6.0.98.gtf.gz > Rattus_norvegicus.Rnor_6.0.98.gtf
     cat Rattus_norvegicus.Rnor_6.0.98.gtf CustomAnnot.gtf > Rattus_norvegicus.Rnor_6.0.98.custom.gtf
     
Generate STAR genomic indexes from custom rat genome and annotation using `genomeGenerate` command from STAR version 2.7.0e with options specified in batch script `STARgenomeGen_custom.sh`  

Generate RSEM-STAR genomic indexes from custom rat genome and annotation using `rsem-prepare-reference` command with RSEM version 1.3.1 with options specified in batch script `RSEM_RnorV98customprep.sh`  

Generate Bowtie2 indexes from Rnor6.0.99 (without custom annotation) using the `bowtie2-build` command with bowtie2 version 2.3.4.2 with options specified in batch script `Bowtie_Rnor6v99_build.sh`

## Part II: Identification of differentially expressed genes  in total RNA using STAR alignment with gene feature counts by HT-Seq and comparison with DEseq2
RNA-seq QC check using `FastQC` version 0.11.5 and illumina NextSeq adapters trimmed using `TrimGalore` version 0.4.5 as specified in `Total_RNA_QCandTRIM.sh` with the sample list file `expList_2017.txt`  

Splice-aware refrence-based alignment using STAR version 2.7.0e was implemented with options specified in `STAR_RnorCustom_mapping.sh` combined with the custom annotation genomic indexes generated in Part I.

Gene-level feature counts were made with HTseq version 0.11.2 using parameters specified in `HTSeq_Rnorcustom_v98.sh`.  
*Note that the settings* `--mode union` *and* `--nonqunique all` *were required to accurately quantify the overlapping bicistronic Snurf-Snrpn transcript and the multicopy Snord115 and Snord116 genes of the PWS domain.*  

The DESeq2 version 3.9 R package (R version 3.6) was used for differential expression anlaysis as implemented in R script `Deseq2_STARHTSeq_RScript.R`.  

*Significance cutoffs of padj<0.10 lenient or padj<0.05 stringent were used in the identification of differentially expressed genes to investigate further.*

## Part III: Alternate identification of differentially expressed genes in total RNA using STAR alignment with RSEM and comparison with DEseq2
Starting with custom genome annotation and pre-processed reads from PartI and PartII. `RSEM` v.1.3.1 with `STAR` was used to perform both alignment and feature counts as implemented in `RSEM_RnorV98custom_calcexp`.  

A normalized matrix of gene-level and isoform-level read-counts for the 6 samples was then generated using `RSEM_RnorV98custom_genedatamatrix.sh`.  

The read-count matrix was subsequently imported into R using `tximport` and differential expression quantified by `DESeq2` as written in the script `RSEM_RnorV98custom_Deseq2`.  

*A Significance cutoff of padj<0.10 was used to compare with results from PartII*


## Part IV: smallRNA-Seq differential expression using Bowtie2 alignment, HTseq and DEseq2

## Part V References and links to software used herein

### STAR
https://github.com/alexdobin/STAR  
https://pubmed.ncbi.nlm.nih.gov/23104886/  
A. Dobin, C.A. Davis, F. Schlesinger, J. Drenkow, C. Zaleski, S. Jha, P. Batut, M. Chaisson, T.R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics 29(1) (2013) 15-21.

### HTSeq
https://github.com/htseq/htseq  
https://htseq.readthedocs.io/en/master/  
https://pubmed.ncbi.nlm.nih.gov/25260700/  
Simon Anders, Paul Theodor Pyl, Wolfgang Huber. HTSeq — A Python framework to work with high-throughput sequencing data Bioinformatics (2014).

### RSEM
http://deweylab.github.io/RSEM/  
https://pubmed.ncbi.nlm.nih.gov/21816040/  


### DESeq2
https://github.com/mikelove/DESeq2  
https://bioconductor.org/packages/DESeq2/  
https://pubmed.ncbi.nlm.nih.gov/25516281/  
Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.    

### Trimomatic
https://github.com/timflutre/trimmomatic  
http://www.usadellab.org/cms/index.php?page=trimmomatic  
https://pubmed.ncbi.nlm.nih.gov/24695404/  
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics.  

# INS-1_PWS_PNAS
bash scripts for RNA-Seq analysis associated with Koppes et al. 2021 for analysis of total- and small- RNA-Seq of control and PWS INS-1 cells

## Part I: Custom Genome Generation
Download Rnor6 top level genome sequence and gene set annotation from Ensembl (https://useast.ensembl.org/info/data/ftp/index.html)  

    wget http://ftp.ensembl.org/pub/release-98/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
    wget http://ftp.ensembl.org/pub/release-98/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.98.gtf.gz

Note that as of 10/21/2021 current genome and annotation is now Rnor_6.0.104  

Unzip the genome fasta file and then append the custom transgne sequences from `mINS2tg.fa` and `hINStg.fa`

     gunzip -c Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz > Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
     cat Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa hINStg.fa mINS2tg.fa > Rattus_norvegicus.Rnor_6.0.dna.custom.fa
     gunzip -c Rattus_norvegicus.Rnor_6.0.104.gtf.gz > Rattus_norvegicus.Rnor_6.0.104.gtf
     cat Rattus_norvegicus.Rnor_6.0.104.gtf CustomAnnot.gtf >

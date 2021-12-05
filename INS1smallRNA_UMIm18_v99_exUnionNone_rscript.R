##Script for DESeq2 analysis of Bowtie2-Htseq small RNA-Seq count data with Cutadapt m=18 exon feature HTseq UnionNone Rnorv99 Annotation
##Erik Koppes 111520
##import and merge all .tsv count data files
library(readr)
library(dplyr)
Con1 <- read_tsv("INS1_line_5_9_plus_S1.tsv", col_names = c("gene_ID", "gene_name","INS1_line_5_9_plus_S1"), comment = "__") %>% select(-2)
Con2 <- read_tsv("INS1_line_2_S2.tsv", col_names = c("gene_ID", "gene_name", "INS1_line_2_S2"), comment = "__") %>% select(-2)
Con3 <- read_tsv("INS1_line_16_S3.tsv", col_names = c("gene_ID", "gene_name","INS1_line_16_S3"), comment = "__") %>% select(-2)
PWS1 <- read_tsv("INS1_line_3_S4.tsv", col_names = c("gene_ID", "gene_name","INS1_line_3_S4"), comment = "__") %>% select(-2)
PWS2 <- read_tsv("INS1_line_19_1_S5.tsv", col_names = c("gene_ID", "gene_name","INS1_line_19_1_S5"), comment = "__") %>% select(-2)
PWS3 <- read_tsv("INS1_line_19_4_S6.tsv", col_names = c("gene_ID", "gene_name","INS1_line_19_4_S6"), comment = "__") %>% select(-2)
EnsemblIDtoName <- read_tsv("INS1_line_5_9_plus_S1.tsv", col_names = c("gene_ID", "gene_name","INS1_line_5_9_plus_S1"), comment = "__") %>% select(-3)

library(tibble)
countData_INS1smallRNA <- Con1 %>% 
  left_join(Con2, by = "gene_ID") %>%
  left_join(Con3, by = "gene_ID") %>%
  left_join(PWS1, by = "gene_ID") %>%
  left_join(PWS2, by = "gene_ID") %>%
  left_join(PWS3, by = "gene_ID") %>%
  column_to_rownames (var = "gene_ID")

summary(countData_INS1smallRNA)

library(ggplot2)
Con1 %>%
  mutate(ADJ = 1 + INS1_line_5_9_plus_S1) %>%
  ggplot () +
  geom_histogram(aes(x = ADJ), stat = "bin", bins = 100) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Raw expression counts") +
  ylab("Number of genes")


##make meta data table
genotype <- c("Con", "Con", "Con", "PWS", "PWS", "PWS")
sampleNames <- c("INS1_line_5_9_plus_S1", "INS1_line_2_S2", "INS1_line_16_S3",
                 "INS1_line_3_S4", "INS1_line_19_1_S5", "INS1_line_19_4_S6")
metaData <- data.frame(genotype, row.names = sampleNames)

all(rownames(metaData) == colnames(countData_INS1smallRNA))

##generate DESeq1 dds object
library(DESeq2)
dds_INS1smallRNA_m18_exUnionNone <- DESeqDataSetFromMatrix(countData = countData_INS1smallRNA,
                                                          colData = metaData,
                                                          design = ~genotype)

##Estimate library size factors
dds_INS1smallRNA_m18_exUnionNone <- estimateSizeFactors(dds_INS1smallRNA_m18_exUnionNone)
sizeFactors(dds_INS1smallRNA_m18_exUnionNone)
sizeFactortbl <- data.frame(sizeFactors(dds_INS1smallRNA_m18_exUnionNone))
write_excel_csv(sizeFactortbl, "INS1smallRNA_UMIm18_v99_exUnionNone_Sizefactors.csv")

normalized_INS1smallRNA_m18_exUnionNone_counts <- data.frame(counts(dds_INS1smallRNA_m18_exUnionNone, normalized = T))
write_excel_csv(normalized_INS1smallRNA_m18_exUnionNone_counts, "INS1smallRNA_m18_v99_exUnionNone_Normcounts.csv")

##sample correlation Plots
vsd_INS1smallRNA_m18_exUnionNone <- vst(dds_INS1smallRNA_m18_exUnionNone, blind = T)
vsd_mat_INS1smallRNA_m18_exUnionNone <- assay(vsd_INS1smallRNA_m18_exUnionNone)
vsd_cor_INS1smallRNA_m18_exUnionNone <- cor(vsd_mat_INS1smallRNA_m18_exUnionNone)
View(vsd_cor_INS1smallRNA_m18_exUnionNone)


library(pheatmap)
pheatmap(vsd_cor_INS1smallRNA_m18_exUnionNone, annotation = metaData)

#sample PCA plot
plotPCA(vsd_INS1smallRNA_m18_exUnionNone, intgroup = "genotype")

##run DESeq
dds_INS1smallRNA_m18_exUnionNone <- DESeq(dds_INS1smallRNA_m18_exUnionNone)

#calc mean ea row
mean_counts <- apply(countData_INS1smallRNA[, 1:6], 1, mean)
variance_counts <- apply(countData_INS1smallRNA[, 1:6], 1, var)
meanvardf <- data.frame (mean_counts, variance_counts)
ggplot(meanvardf) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene")

##plot dispersion estimates
plotDispEsts(dds_INS1smallRNA_m18_exUnionNone)

##PlotMA
results(dds_INS1smallRNA_m18_exUnionNone, alpha = 0.05)
INS1smallRNA_m18_exUnionNone_res <- results(dds_INS1smallRNA_m18_exUnionNone,
                                           contrast = c("genotype", "PWS", "Con"),
                                           alpha = 0.05)
plotMA(INS1smallRNA_m18_exUnionNone_res, ylim = c(-10, 10))

INS1smallRNA_m18_exUnionNone_res_lfc <- lfcShrink(dds_INS1smallRNA_m18_exUnionNone,
                                                  coef = 2,
                                                  type = "apeglm",)
plotMA(INS1smallRNA_m18_exUnionNone_res_lfc, ylim = c(-5, 5))


INS1smallRNA_m18_exUnionNone_PWS_res_all <- data.frame(INS1smallRNA_m18_exUnionNone_res_lfc) %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(EnsemblIDtoName, by = "gene_ID") %>%
  select("gene_ID", "gene_name", everything())

INS1smallRNA_m18_exUnionNone_PWS_res_all_wnormcounts <- normalized_INS1smallRNA_m18_exUnionNone_counts  %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(INS1smallRNA_m18_exUnionNone_PWS_res_all) %>%
  select("gene_ID", "gene_name", everything())

INS1smallRNA_m18_exUnionNone_PWS_sig0.05 <- INS1smallRNA_m18_exUnionNone_PWS_res_all_wnormcounts %>%
  filter(padj <= 0.05) %>%
  arrange(padj)

INS1smallRNA_m18_exUnionNone_PWS_sig0.1 <- INS1smallRNA_m18_exUnionNone_PWS_res_all_wnormcounts %>%
  filter(padj <= 0.10) %>%
  arrange(padj)

write_excel_csv(data.frame(INS1smallRNA_m18_exUnionNone_PWS_res_all_wnormcounts), "INS1smallRNA_UMIm18_v99_exUnionAll_PWS_res_all_wnormcounts.csv")
write_excel_csv(data.frame(INS1smallRNA_m18_exUnionNone_PWS_sig0.05), "INS1smallRNA_UMIm18_exUnionAll_v99_PWS_sig0.05.csv")
write_excel_csv(data.frame(INS1smallRNA_m18_exUnionNone_PWS_sig0.1), "INS1smallRNA_UMIm18_exUnionAll_v99_PWS_sig0.1.csv")
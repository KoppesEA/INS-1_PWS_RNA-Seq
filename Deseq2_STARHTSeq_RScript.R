##Script for DESeq2 analysis of STAR-Htseq (Mode:Union + Nonunique:all) count data with custom Rnor v98 annotation
##Erik Koppes 030520
##import and merge all .tsv count data files
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

Con1 <- read_tsv("Control_1_S1.tsv", col_names = c("gene_ID", "gene_name","Control_1_S1"), comment = "__") %>% select(-2)
Con2 <- read_tsv("Control_2_S2.tsv", col_names = c("gene_ID", "gene_name","Control_2_S2"), comment = "__") %>% select(-2)
Con3 <- read_tsv("Control_3_S3.tsv", col_names = c("gene_ID", "gene_name","Control_3_S3"), comment = "__") %>% select(-2)
PWS1 <- read_tsv("PWS_1_S4.tsv", col_names = c("gene_ID", "gene_name","PWS_1_S4"), comment = "__") %>% select(-2)
PWS2 <- read_tsv("PWS_2_S5.tsv", col_names = c("gene_ID", "gene_name","PWS_2_S5"), comment = "__") %>% select(-2)
PWS3 <- read_tsv("PWS_3_S6.tsv", col_names = c("gene_ID", "gene_name","PWS_3_S6"), comment = "__") %>% select(-2)
EnsemblIDtoName <- read_tsv("Control_1_S1.tsv", col_names = c("gene_ID", "gene_name","Control_1_S1"), comment = "__") %>% select(-3)

countData_INS1custom <- Con1 %>% 
  left_join(Con2, by = "gene_ID") %>%
  left_join(Con3, by = "gene_ID") %>%
  left_join(PWS1, by = "gene_ID") %>%
  left_join(PWS2, by = "gene_ID") %>%
  left_join(PWS3, by = "gene_ID") %>%
  column_to_rownames (var = "gene_ID")

Con1 %>%
  mutate(ADJ = 1 + Control_1_S1) %>%
  ggplot () +
  geom_histogram(aes(x = ADJ), stat = "bin", bins = 100) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Raw expression counts") +
  ylab("Number of genes")

##make meta data table
genotype <- c("Con", "Con", "Con", "PWS", "PWS", "PWS")
sampleNames <- c("Control_1_S1", "Control_2_S2", "Control_3_S3", "PWS_1_S4", "PWS_2_S5", "PWS_3_S6")
metaData <- data.frame(genotype, row.names = sampleNames)

all(rownames(metaData) == colnames(countData_INS1custom))

library(DESeq2)
dds_INS1custom <- DESeqDataSetFromMatrix(countData = countData_INS1custom,
                                         colData = metaData,
                                         design = ~genotype)

dds_INS1custom <- estimateSizeFactors(dds_INS1custom)
sizeFactors(dds_INS1custom)
sizeFactortbl <- data.frame(sizeFactors(dds_INS1custom))
write_excel_csv(sizeFactortbl, "INS1cutstomSizefactors_UnionAll.csv")

normalized_INS1custom_counts <- data.frame(counts(dds_INS1custom, normalized = T))
write_excel_csv(normalized_INS1custom_counts, "INS1cutstomNormcounts_UnionAll.csv")

##sample correlation Plots
vsd_INS1 <- vst(dds_INS1custom, blind = T)
vsd_mat_INS1 <- assay(vsd_INS1)
vsd_cor_INS1 <- cor(vsd_mat_INS1)
View(vsd_cor_INS1)

library(pheatmap)
pheatmap(vsd_cor_INS1, annotation = metaData)

#sample PCA plot
plotPCA(vsd_INS1, intgroup = "genotype")

##run DESeq
dds_INS1custom <- DESeq(dds_INS1custom)

#calc mean ea row
mean_counts <- apply(countData_INS1custom[, 1:6], 1, mean)
variance_counts <- apply(countData_INS1custom[,1:6], 1, var)
meanvardf <- data.frame (mean_counts, variance_counts)
ggplot(meanvardf) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene")

##plot dispersion estimates
plotDispEsts(dds_INS1custom)

results(dds_INS1custom, alpha = 0.05)
INS1_PWS_res <- results(dds_INS1custom,
                        contrast = c("genotype", "PWS", "Con"),
                        alpha = 0.05)
plotMA(INS1_PWS_res, ylim = c(-10, 10))

INS1_PWS_res_all <- data.frame(INS1_PWS_res) %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(EnsemblIDtoName, by = "gene_ID") %>%
  select("gene_ID", "gene_name", everything())

INS1_PWS_res_all_wnormcounts <- normalized_INS1custom_counts  %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(INS1_PWS_res_all) %>%
  select("gene_ID", "gene_name", everything())

INS1_PWS_sig0.05 <- INS1_PWS_res_all_wnormcounts %>%
  filter(padj <= 0.05) %>%
  arrange(padj)

INS1_PWS_sig0.1 <- INS1_PWS_res_all_wnormcounts %>%
  filter(padj <= 0.10) %>%
  arrange(padj)

write_excel_csv(data.frame(INS1_PWS_res_all_wnormcounts), "INS1_PWS_res_all_v98cust_UnionAll.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.05), "INS1_PWS_res_padj_v98cust_UnionAll.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.1), "INS1_PWS_res_padj0.1_v98cust_UnionAll.csv")

##Analysis with lfc_Shrink 
##Don't use as final as shrinks logfold change of PWS genes too much
INS1_PWS_res_lfcShrink <- lfcShrink(dds_INS1custom,
                                    contrast = c("genotype", "PWS", "Con"),
                                    alpha = 0.05)
plotMA(INS1_PWS_res_lfcShrink, ylim = c(-5, 5))

INS1_PWS_res_all_lfcShrink <- data.frame(INS1_PWS_res_lfcShrink) %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(EnsemblIDtoName, by = "gene_ID") %>%
  select("gene_ID", "gene_name", everything())

INS1_PWS_res_all_wnormcounts_lfcShrink <- normalized_INS1custom_counts  %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(INS1_PWS_res_all) %>%
  select("gene_ID", "gene_name", everything())

INS1_PWS_sig0.05_lfcShrink <- INS1_PWS_res_all_wnormcounts_lfcShrink %>%
  filter(padj <= 0.05) %>%
  arrange(padj)

INS1_PWS_sig0.1_lfcShrink <- INS1_PWS_res_all_wnormcounts_lfcShrink %>%
  filter(padj <= 0.10) %>%
  arrange(padj)

write_excel_csv(data.frame(INS1_PWS_res_all_wnormcounts_lfcShrink), "INS1_PWS_res_all_v98cust_UnionAll_lfcShrink.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.05_lfcShrink), "INS1_PWS_res_padj_v98cust_UnionAll_lfcShrink.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.1_lfcShrink), "INS1_PWS_res_padj0.1_v98cust_UnionAll_lfcShrink.csv")

##Draw Volcano Plot
INS1_PWS_res_all %>%
  mutate(threshold = padj <= 0.05) %>%
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 2)) +
  scale_x_continuous(limits = c(-12.0, 6.0), breaks = seq(-12.0, 6.0, 3.0), expand = c(0,0)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(linetype = 2, aes(xintercept = 0))


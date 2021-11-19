##Script for DESeq2 analysis of STAR-RSEM count data with custom Rnor v98 annotation
##Erik Koppes 040720
##import and merge all .tsv count data files


library(tximportData)
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

INS1_GeneMat <- data.matrix(read.table(file="INS1_RSEM_genes.matrix"))
colnames(INS1_GeneMat) <- c("Con_1", "Con_2", "Con_3", "PWS_1", "PWS_2", "PWS_3")

dir_results <- "/bgfs/rnicholls/hpark_share/EK_100618/RSEM_Rnor_Custom_V98"
list.files(dir_results)
sample <- c("Control_1_S1", "Control_2_S2", "Control_3_S3",
            "PWS_1_S4", "PWS_2_S5", "PWS_3_S6")
genotype <- c("Con", "Con", "Con", "PWS", "PWS", "PWS")
samples <- data.frame(cbind(sample, genotype))

##from vignette; strings asfactor import
#dir <- system.file("extdata", package = "tximportData")
#samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

##in vignente it is samples$run but here run=sample since no duplicates
files <- file.path(dir_results, paste0(samples$sample, ".genes.results"))
names(files) <-  c("Con_1", "Con_2", "Con_3", "PWS_1", "PWS_2", "PWS_3")
txi_rsem <- tximport(files, type ="rsem", txIn = FALSE, txOut = FALSE)

#after running DESeqDataSetFromTximport need to set any of txi_rsem$length that are zero to equal 1
txi_rsem$length[txi_rsem$length == 0] <- 1

sampleTable <- data.frame(genotype = factor(rep(c("Con", "PWS"), each = 3)))
rownames(sampleTable) <- colnames(txi_rsem$counts)

dds_INS1rsem <- DESeqDataSetFromTximport(txi_rsem, sampleTable, ~genotype)

dds_INS1rsem <- estimateSizeFactors(dds_INS1rsem)
sizeFactors(dds_INS1rsem)
sizeFactortbl <- data.frame(sizeFactors(dds_INS1rsem))
write_excel_csv(sizeFactortbl, "INS1cutstomSizefactors_RSEM.csv")

normalized_INS1custom_counts <- data.frame(counts(dds_INS1rsem, normalized = T))
write_excel_csv(normalized_INS1custom_counts, "INS1cutstomNormcounts_RSEM.csv")

##sample correlation Plots
vsd_INS1 <- vst(dds_INS1rsem, blind = T)
vsd_mat_INS1 <- assay(vsd_INS1)
vsd_cor_INS1 <- cor(vsd_mat_INS1)
View(vsd_cor_INS1)

library(pheatmap)
pheatmap(vsd_cor_INS1, annotation = sampleTable)

#sample PCA plot
plotPCA(vsd_INS1, intgroup = "genotype")

##run DESeq
dds_INS1rsem <- DESeq(dds_INS1rsem)

#calc mean ea row
countData_INS1_RSEM <- txi_rsem$counts
mean_counts <- apply(countData_INS1_RSEM[, 1:6], 1, mean)
variance_counts <- apply(countData_INS1_RSEM[,1:6], 1, var)
meanvardf <- data.frame (mean_counts, variance_counts)
ggplot(meanvardf) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene")

##plot dispersion estimates
plotDispEsts(dds_INS1rsem)

results(dds_INS1rsem, alpha = 0.05)
INS1_PWS_res <- results(dds_INS1rsem,
                        contrast = c("genotype", "PWS", "Con"),
                        alpha = 0.05)
plotMA(INS1_PWS_res, ylim = c(-10, 10))

INS1_PWS_res_all <- data.frame(INS1_PWS_res) %>%
  rownames_to_column(var = "gene_ID") %>%
  right_join(Ensembl_IDtoGene, by = "gene_ID") %>%
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

write_excel_csv(data.frame(INS1_PWS_res_all_wnormcounts), "INS1_PWS_res_all_v98cust_RSEM_DESEQ2.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.05), "INS1_PWS_res_padj_v98cust_RSEM_DESEQ2.csv")
write_excel_csv(data.frame(INS1_PWS_sig0.1), "INS1_PWS_res_padj0.1_v98cust_RSEM_DESEQ2.csv")

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


#preprocessing
library(readr)
library(dplyr)
library(ggpubr)
source("tools.R")
source("signature_lib.R")
#expression
exp_matrix = read_tsv("./TCGA-KIRC.htseq_fpkm.tsv", col_names = T)
exp_matrix = data.frame(exp_matrix, check.names = F)
#gene annotation
gene_anno = read_tsv("gencode.v22.annotation.gene.probeMap", col_names = T)
gene_anno = data.frame(gene_anno, check.names = F)
gene_anno = gene_anno[,c("id", "gene")]
#join gene anno with expression matrix
exp_matrix_gene_anno = inner_join(gene_anno, exp_matrix, by=c("id"="Ensembl_ID"))
exp_matrix_gene_anno = exp_matrix_gene_anno[!duplicated(exp_matrix_gene_anno[c("gene")]),]
rownames(exp_matrix_gene_anno) = exp_matrix_gene_anno$gene
exp_matrix_gene_anno = exp_matrix_gene_anno[,3:ncol(exp_matrix_gene_anno)]
clinical_df = read_tsv("TCGA-KIRC.GDC_phenotype.tsv", col_names = T)
clinical_df = data.frame(clinical_df, check.names = F)[,c("submitter_id.samples","sample_type.samples")]
clinical_df = clinical_df[clinical_df$sample_type.samples=="Primary Tumor",]


#preprocessing
library(readr)
library(dplyr)
source("tools.R")
#expression
prepare_gene_groups = function(gene_name){
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
gene_name = c("MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8",
              "MT-ND1","MT-ND3","MT-ND4","MT-ND4L")
if(all(gene_name %in% rownames(exp_matrix_gene_anno))){
  exp_vector = exp_matrix_gene_anno[gene_name, ]
}else{
  print(paste0(gene_name, " is not available for analysis"))
}
exp_vector["MT", ] = apply(exp_vector, MARGIN = 2, FUN = mean)
exp_vector = exp_vector[nrow(exp_vector),]

exp_vector_t = data.frame(t(exp_vector), check.names = F)
exp_vector_t$sample_id = rownames(exp_vector_t)
clinical_df = read_tsv("TCGA-KIRC.GDC_phenotype.tsv", col_names = T)
clinical_df = data.frame(clinical_df, check.names = F)[,c("submitter_id.samples","sample_type.samples")]
clinical_df = clinical_df[clinical_df$sample_type.samples=="Primary Tumor",]
exp_vector_t = exp_vector_t[exp_vector_t$sample_id %in% clinical_df$submitter_id.samples,]

survival_raw = data.frame(read_tsv("TCGA-KIRC.survival.tsv", col_names = T),check.names = F)
survival_data = inner_join(exp_vector_t, survival_raw, by=c("sample_id"="sample"))
survival_data = survival_data[!duplicated(survival_data[c("_PATIENT")]),]
survival_data = survival_data %>% mutate(OS.time = round(OS.time/30.417, digit=0))
marker = "MT"
threshold = cutoff_function_options(survival_data, marker, "mean",
                                    custom = 1.5, "OS.time", "OS")
group_var = paste0("gene_expression")
survival_data[[group_var]] <- sapply(survival_data[[marker]], function(x){
  return(ifelse(x>=threshold,"High","Low"))
})

png("survival_test.png", units = "in", res = 300, width = 12, height = 9)
run_survival(survival_data, group_var, marker, threshold)
dev.off()





}




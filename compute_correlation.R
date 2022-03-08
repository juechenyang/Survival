library(readr)
library(dplyr)
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

gene_list = list("MT" = mt_gene
                 ,"Complex_I"=Complex_I
                 ,"Complex_II"=Complex_II
                 ,"Complex_III"=Complex_III
                 ,"Complex_IV"=Complex_IV
                 ,"Complex_V"=Complex_V
                 ,"glycolyticAndTCA"=glycolyticAndTCA
                 ,"MRP_positive"=MRP_positive
                 ,"MRP_negative"=MRP_negative)

x=9
gene_signature_name = names(gene_list)[x]
gene_names = unlist(gene_list[gene_signature_name], use.names = F)
is_in_exp = gene_names %in% rownames(exp_matrix_gene_anno)
if(!all(is_in_exp)){
  not_available_genes = gene_names[!is_in_exp]
  print(paste0(paste(not_available_genes, collapse = ","), 
               " is not available for analysis"))
}
gene_names = gene_names[is_in_exp]
exp_vector = exp_matrix_gene_anno[gene_names, ]
exp_vector[gene_signature_name, ] = apply(exp_vector, MARGIN = 2, FUN = mean)
exp_vector = exp_vector[nrow(exp_vector),]
exp_vector_t = data.frame(t(exp_vector), check.names = F)
exp_vector_t$sample_id = rownames(exp_vector_t)
exp_vector_t = exp_vector_t[exp_vector_t$sample_id %in% clinical_df$submitter_id.samples,]

selected_samples = names(exp_matrix_gene_anno)[names(exp_matrix_gene_anno) %in% exp_vector_t$sample_id]
exp_matrix_gene_anno_bk = exp_matrix_gene_anno[,selected_samples]
correlation_value = apply(exp_matrix_gene_anno_bk, MARGIN = 1, FUN=function(x){
  cor_model = cor.test(as.numeric(x), exp_vector_t[,gene_signature_name],
                       method = "spearman")
  return(c(cor_model$estimate, cor_model$p.value))
})
correlation_value = data.frame(t(correlation_value))
names(correlation_value) = c("cor", "p.value")
correlation_value = na.omit(correlation_value)
correlation_value = correlation_value[order(correlation_value[,"cor"], decreasing = T), ]
write.csv(correlation_value, paste0(gene_signature_name, "_correlated_genes.csv"))


#compute 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
ids<-bitr(rownames(correlation_value), fromType = "SYMBOL", 
          toType = "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
correlation_value = correlation_value[rownames(correlation_value) %in% dedup_ids$SYMBOL,]
kegg_gene_list = correlation_value$cor
names(kegg_gene_list) = dedup_ids$ENTREZID
kegg_gene_list = sort(kegg_gene_list, decreasing = T)

GSEA_kegg <- gseKEGG(geneList = kegg_gene_list,
                     organism     = "hsa",
                     minGSSize    = 10,
                     maxGSSize    = 1000,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "none",
                     keyType       = "ncbi-geneid")
png(paste0(gene_signature_name, "_gsea_dot_plot.png"), units = "in", res = 300, width = 8, height = 12)
dotplot(GSEA_kegg, showCategory = 10, 
        title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()
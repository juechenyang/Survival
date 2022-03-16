#preprocessing
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
clinical_features = c("submitter_id.samples","sample_type.samples","gender.demographic")
clinical_df = data.frame(clinical_df, check.names = F)[,clinical_features]
clinical_df = clinical_df[clinical_df$sample_type.samples=="Primary Tumor",]
#clinical_df = clinical_df[clinical_df$gender.demographic=="male",]
survival_raw = data.frame(read_tsv("TCGA-KIRC.survival.tsv", col_names = T),check.names = F)

get_survival_plot = function(gene_list, cohort_25 = F, percent=0.25){
  plot_list = list()
  for(x in 1:length(gene_list)){
    gene_signature_name = names(gene_list)[x]
    gene_names = unlist(gene_list[gene_signature_name], use.names = F)
    is_in_exp = gene_names %in% rownames(exp_matrix_gene_anno)
    if(FALSE %in% is_in_exp){
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
    
    
    survival_data = inner_join(exp_vector_t, survival_raw, by=c("sample_id"="sample"))
    survival_data = survival_data[!duplicated(survival_data[c("_PATIENT")]),]
    survival_data = survival_data %>% mutate(OS.time = round(OS.time/30.417, digit=0))
    #keep data with OS less than 10 years
    #survival_data = survival_data[survival_data$OS.time<=120,]
    threshold = cutoff_function_options(survival_data, gene_signature_name, "mean",
                                        custom = 1.5, "OS.time", "OS")
    group_var = paste0(gene_signature_name, "_group")
    if(cohort_25){
      top25_threshold = quantile(survival_data[[gene_signature_name]], (1-percent))
      bottom25_threshold = quantile(survival_data[[gene_signature_name]], percent)
      survival_data[[group_var]] <- sapply(survival_data[[gene_signature_name]], function(x){
        if(x >= top25_threshold){
          return("top25")
        }else if(x <= bottom25_threshold){
          return("bottom25")
        }
        else{
          return("exclude")
        }
      })
      survival_data = survival_data[!survival_data[[group_var]]=="exclude",]
      threshold = c(top25_threshold, bottom25_threshold)
    }else{
      survival_data[[group_var]] <- sapply(survival_data[[gene_signature_name]], function(x){
        return(ifelse(x>=threshold,"High","Low"))
      })
    }
    survival_plot = run_survival(survival_data, group_var, gene_signature_name, threshold)
    plot_list[[x]] = survival_plot
  }
  return(plot_list)
}


gene_list = list("MT" = MT
                 ,"Complex_I"=Complex_I
                 ,"Complex_II"=Complex_II
                 ,"Complex_III"=Complex_III
                 ,"Complex_IV"=Complex_IV
                 ,"Complex_V"=Complex_V
                 ,"TCA"=TCA
                 ,"glycolysis"=glycolysis
#                 ,"glycolyticAndTCA"=glycolyticAndTCA
                 ,"MAS"=MAS
                 ,"HIF1A" = HIF1A
#                 ,"CU"=CU
                 ,"EMT"=EMT
#                 ,"HIF1A_2A"=HIF1A_2A
                 ,"HLA"=HLA
                 ,"MRP"=MRP
                 ,"MRP_positive"=MRP_positive
                 ,"MRP_negative"=MRP_negative
                 )
mat_seq = matrix(names(gene_list),nrow = 5)
mat_seq = t(mat_seq)
mat_seq_v = as.character(mat_seq)
gene_list = gene_list[mat_seq_v]
plot_list = get_survival_plot(gene_list = gene_list, cohort_25 = F)




png("survival_mean_compare.png", units = "in", res = 400, width = 20, height = 12)
arrange_ggsurvplots(plot_list, print = TRUE,
                    ncol = 5, nrow = 3
#                    ,risk.table.height = 0.3
                    )
dev.off()

png("survival_25percent_compare.png", units = "in", res = 400, width = 18, height = 11)
arrange_ggsurvplots(plot_list, print = TRUE,
                    ncol = 5, nrow = 3
#                   ,risk.table.height = 0.3
                    )
dev.off()

# gene_list = list("new_sig"=new_signature)
# plot_list = get_survival_plot(gene_list = gene_list, cohort_25 = T)
# plt_list[[2]] =  plot_list[[1]]
# png("survival_new_sig.png", units = "in", res = 300, width = 10, height = 6)
# arrange_ggsurvplots(plt_list, print = TRUE,
#                     ncol = 2, risk.table.height = 0.3)
# dev.off()



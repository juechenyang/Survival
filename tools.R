library(survival)
library(survminer)
library(maxstat)
library(ggplot2)
cutoff_function_default = function(dataframe,input, time, status){
  res.cut = surv_cutpoint(dataframe,time=time,event=status,variables=input)
  return(res.cut)
}
#Different cutoff options
cutoff_function_options = function(dataframe,input,cutoption,custom=NULL, time=NULL, status=NULL){
  if(cutoption == "auto_calculate"){
    cutvalue = summary(cutoff_function_default(dataframe,input, time, status))[1,1]
    return(cutvalue)
  }else if(cutoption == "median"){
    cutvalue = median(dataframe[[input]])
  }else if(cutoption == "mean"){
    cutvalue = mean(dataframe[[input]])
  }else if(cutoption == "custom"){
    cutvalue = custom
  }
  return(cutvalue)
}

get_surv_obj <- function(df){
  surv_object = Surv(time=df$OS.time,event=df$OS)
  return(surv_object)
}

run_survival <- function(df, group, parameter, threshold){
  surv_object <- get_surv_obj(df)
  fit_geneexpr = surv_fit(as.formula(paste('surv_object ~', group)),data=df)
  group_tags = unique(df[[group]])
  if("top25" %in% group_tags){
    colors = c("top25"="#f9837b", "bottom25"="#09c1c6")
    palette_colors = as.character(colors[c("bottom25", "top25")])
    group_text = c("bottom25", "top25")
  }else{
    colors = c("High"="#f9837b", "Low"="#09c1c6")
    palette_colors = as.character(colors[c("High", "Low")])
    group_text = c("High", "Low")
  }
  
  if(grepl("_", parameter)){
    parameter = stringr::str_replace_all(parameter, "_", " ")
  }
  plot_fit = ggsurvplot(fit_geneexpr,data=df
#                        ,tables.theme = theme_survminer(font.tickslab = 10)
                        ,risk.table = F
#                        ,risk.table.height=0.44
                        ,fontsize = 50,xlab = "Months"
                        ,legend = "none"
                        ,legend.lab = group_text
                        ,legend.title = parameter
                        ,palette = palette_colors
                        ,xlim=c(0,120)
                        )
  fit_pmodel = surv_pvalue(fit_geneexpr,data=df,method="survdiff")
  fit_pval=fit_pmodel$pval
  
  # if("High" %in% group_tags){
  #   high_count = table(df[[group]])[group_tags[1]]
  #   low_count = table(df[[group]])[group_tags[2]]
  #   fit_annot = sprintf("high-low groups cutoff = %.2f\nP-value(Log-Rank) = %.5f\nHigh group count = %s\nLow group count = %s",threshold,fit_pval,high_count, low_count)
  # }else{
  #   high_count = table(df[[group]])[group_tags[1]]
  #   low_count = table(df[[group]])[group_tags[2]]
  #   threshold = round(threshold, 2)
  #   threshold = paste(threshold, collapse = ",")
  #   fit_annot = sprintf("groups cutoff = %s\nP-value(Log-Rank) = %.5f\ntop25 group count = %s\nbottom group count = %s",threshold,fit_pval,high_count, low_count)
  # }
  fit_annot = sprintf("P = %.5f",fit_pval)
  plot_fit$plot = plot_fit$plot+
    annotate("text",x=0,y=0.3,label=fit_annot,size=5,hjust=0)+
    scale_x_continuous(breaks = seq(0,120,30))+
    ggtitle(parameter)+
    theme(plot.title = element_text(size=15, face="bold", hjust=0.5))
  return(plot_fit)
}
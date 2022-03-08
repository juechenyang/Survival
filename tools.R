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
  }else{
    colors = c("High"="#f9837b", "Low"="#09c1c6")
    palette_colors = as.character(colors[c("High", "Low")])
  }
  group_text = paste0(parameter, " ", group_tags)
  plot_fit = ggsurvplot(fit_geneexpr,data=df,risk.table = TRUE,
                        tables.theme = theme_survminer(font.tickslab = 10), 
                        risk.table.height=0.44, fontsize = 3,xlab = "Time in Months"
#                        ,legend.lab = group_text
                        ,palette = palette_colors
                        )
  fit_pmodel = surv_pvalue(fit_geneexpr,data=df,method="survdiff")
  fit_pval=fit_pmodel$pval
  
  if("High" %in% group_tags){
    high_count = table(df[[group]])[group_tags[1]]
    low_count = table(df[[group]])[group_tags[2]]
    fit_annot = sprintf("high-low groups cutoff = %.2f\nP-value(Log-Rank) = %.5f\nHigh group count = %s\nLow group count = %s",threshold,fit_pval,high_count, low_count)
  }else{
    high_count = table(df[[group]])[group_tags[1]]
    low_count = table(df[[group]])[group_tags[2]]
    threshold = round(threshold, 2)
    threshold = paste(threshold, collapse = ",")
    fit_annot = sprintf("groups cutoff = %s\nP-value(Log-Rank) = %.5f\ntop25 group count = %s\nbottom group count = %s",threshold,fit_pval,high_count, low_count)
  }
  plot_fit$plot = plot_fit$plot+annotate("text",x=0,y=0.3,label=fit_annot,size=3,hjust=0)
  return(plot_fit)
}
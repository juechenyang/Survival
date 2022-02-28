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

run_survival <- function(df, group, parameter, threshold=NULL){
  surv_object <- get_surv_obj(df)
  fit_geneexpr = surv_fit(as.formula(paste('surv_object ~', group)),data=df)
  group_text = c(paste0(parameter, " high"), paste0(parameter, " low"))
  plot_fit = ggsurvplot(fit_geneexpr,data=df,risk.table = TRUE,
                        tables.theme = theme_survminer(font.tickslab = 10), 
                        risk.table.height=0.44, fontsize = 3,xlab = "Time in Months",
                        legend.lab = group_text)
  fit_pmodel = surv_pvalue(fit_geneexpr,data=df,method="survdiff")
  fit_pval=fit_pmodel$pval
  
  if(class(df[[parameter]])=='numeric'){
    high_count = table(df[[group]])[c('High')]
    low_count = table(df[[group]])[c('Low')]
    fit_annot = sprintf("high-low groups cutoff = %.2f\nP-value(Log-Rank) = %.5f\nHigh group count = %s\nLow group count = %s",threshold,fit_pval,high_count, low_count)
  }else{
    fit_annot = sprintf("P-value(Log-Rank) = %.5f",fit_pval)
    all_groups = table(df[[group]])
    anno_str=''
    for(i in 1:length(all_groups)){
      anno_str = paste0(anno_str, '#', names(all_groups)[[i]], '=', as.character(all_groups[[i]]), '\n')
    }
    count_annot=anno_str
  }
  plot_fit$plot = plot_fit$plot+annotate("text",x=0,y=0.3,label=fit_annot,size=3,hjust=0)
  return(plot_fit)
}
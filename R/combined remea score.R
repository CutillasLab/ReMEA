
combined_remea_score <- function(df.remea_results){


  .get_mean_and_pvalue <- function(dfr,by="delta_rank_median"){

    xx <- reshape2::dcast(perturbagens~ signature.type,value.var =  by,
                             data=dfr)
    nums <- unlist(lapply(xx, is.numeric), use.names = FALSE)
    #xx$pval <- apply( scale(xx[,nums]),1, function(x) if(length(x)>3){ t.test(x,alternative="greater",mu=0)$p.value})
    xx$mean <- scale(apply( xx[,nums],1,mean,na.rm=T))
    #xx$qvalue <- p.adjust(xx$pval,method = "fdr")
    #xx$alpha <- (-log10(xx$pval)) * xx$mean
    return(xx)
  }

  dfx.med <- .get_mean_and_pvalue(dfr=df.remea_results,"delta_rank_median")
  dfx.mean<- .get_mean_and_pvalue(dfr=df.remea_results,"delta_rank_mean")
  dfx.geommean<- .get_mean_and_pvalue(dfr=df.remea_results,"delta_rank_geommean")
  dfx.zscore<- .get_mean_and_pvalue(dfr=df.remea_results,"delta.zscore")

  df.combi <- data.frame(perturbagen=dfx.zscore$perturbagens,
                         #av.p = apply(data.frame(dfx.zscore$qvalue,dfx.geommean$qvalue,dfx.mean$qvalue,dfx.med$qvalue),1,min),
                          av.effect = apply(data.frame(dfx.zscore$mean,dfx.geommean$mean,dfx.mean$mean,dfx.med$mean),1,max))


  df.combi$combined_pvalue <- ReMEA::combined_pvalues(df.remea_results)
  df.combi$alpha <- -log10(df.combi$combined_pvalue)* df.combi$av.effect
  df.combi$qvalue <- p.adjust(df.combi$combined_pvalue,method = "fdr")



  sig.combi <- subset(df.combi,df.combi$qvalue<0.05 & df.combi$av.effect>0)
  sign.per <- sig.combi$perturbagen
  length(sign.per)

  x <- sig.combi[order(sig.combi$combined_pvalue),]
 #x <- sig.combi[order(-sig.combi$alpha),]
  #x[x$perturbagen=="PDPK1",]
  #x[x$perturbagen=="PIK3CA",]
  #x[x$perturbagen=="BCL2",]
  #x[x$perturbagen=="MTOR",]
  #x[x$perturbagen=="IRS1",]

  #df.combi[df.combi$perturbagen=="AKT2",]

  pplot <- ggplot(x[1:80,],aes(x=av.effect,y=reorder(perturbagen,av.effect),color=-log10(qvalue),size=-log10(qvalue)))+
    geom_point()+
    scale_color_gradient2(low="orange",mid = "orange",high = "red")+
    theme_bw()

  xx <- 0
  if (xx==1){
    pe <- ReMEA::pathway_enrichment_remea(genes = sign.per,
                                        background.genes = df.combi$perturbagen,
                                        prot_dbs = c("kegg","hallmark.genes","nci","process"))

    pe <- pe[pe$counts>3,]
    pe <- pe[order(pe$pvalues),]
 }

  return(list(combined_score=df.combi,
              significant_genes=sign.per,
              remea_plot=pplot))
  }


pcom_by_stouffer <- function(p){
  .erf <- function(x) 2 * pnorm(2 * x/ sqrt(2)) - 1
  .erfinv <- function(x) qnorm( (x+1)/2 ) / sqrt(2)
  .pcomb <- function(p) (1-.erf(sum(sqrt(2) * .erfinv(1-2*p))/sqrt(2*length(p))))/2
  p <- na.omit(unlist(p))
  pl <- NA
  pl <- length(p)
  { if (is.na(pl)) { res <- "There was an empty array of p-values"}
    else
      res <- .pcomb(p) }
  return(res)
}


combined_pvalues<- function(df.remea_results){
  xx <- reshape2::dcast(perturbagens~ signature.type,value.var = "pvalue.ks",data=df.remea_results)
  xx <- xx[,2:ncol(xx)]
  pvals <- apply(xx, 1,  pcom_by_stouffer)
  #xx$pvals <- pvals
  return(pvals)
}

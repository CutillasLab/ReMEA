

find_top_peturbagens_in_ReMEA_output <- function(df.ReMEA_ouput, top_n=40, counts_cutoff=15){

  df <- df.ReMEA_ouput
  df <- subset(df,df$total_counts>counts_cutoff)
  dfx <- reshape2::dcast(perturbagens~ signature.type,value.var =  'max_delta_score',data=df)
  nums <- unlist(lapply(dfx, is.numeric), use.names = FALSE)

  n_measurements <- ncol(dfx[,nums])

  dfx$count <- apply( dfx[,nums],1, function(x) sum(x>0, na.rm = T))
  dfx$mean <- apply( dfx[,nums],1,mean,na.rm=T)
  dfx <- subset(dfx,dfx$count>round(n_measurements/2))
  dfx <- dfx[order(-dfx$mean),]

  if (nrow(dfx)>top_n){
  top_per <- dfx$perturbagens[1:top_n]
  }else{
    top_per <- dfx$perturbagens
  }
  return(na.omit(top_per))

  }

compare_remea_vs_prot_expression <- function(df.remea, df.protein){

  df.remea <- df
  df.protein <- protein_data

  head(df.remea)
  head(df.protein)

  dfxx <- reshape2::dcast(perturbagens~ signature.type,value.var =  'delta_rank',data=df.remea)
  nums <- unlist(lapply(dfxx, is.numeric), use.names = FALSE)
  dfxx$count <- apply( dfxx[,nums],1, function(x) sum(x>0, na.rm = T))
  dfxx$mean <- apply( dfxx[,nums],1,mean,na.rm=T)
  head(dfxx)

  df.protein$protein <- gsub(";","",df.protein$protein,fixed = T)

  colnames(df.protein)[2] <- "Expression"

  dfx <- merge.data.frame(dfxx, df.protein,by.x = "perturbagens",by.y="protein")

  #dfx <- subset(dfx,dfx$total_counts>1 & dfx$qvalue<0.05 & dfx$FDR<0.05)

  head(dfx)

  sig.prots <- df.protein[df.protein$FDR<0.1 & df.protein$Expression>0,]$protein
  all.prots <- df.protein$protein

  e.prot <- pathway_enrichment_remea(genes=sig.prots,
                                     background.genes = all.prots,
                                     prot_dbs=c("kegg","hallmark.genes","nci","process","selected"))

  sig.pers <- df.remea[df.remea$qvalue<0.05 & df.remea$delta_rank>0,]$perturbagens
  all.pers <- df.remea$perturbagens


  e.pers <- pathway_enrichment_remea(genes=sig.pers,
                                     background.genes = all.pers,
                                     prot_dbs=c("kegg","hallmark.genes","nci","process","selected"))

  head(e.pers)



  ee <- merge.data.frame(e.prot,e.pers,by="pathway")

  ee$combined_pvalue <- ee$FDR.x*ee$FDR.y

  ees <- subset(ee,ee$pvalues.x<0.05 & ee$pvalues.y<0.05)

  ggplot(ees,aes(x=log2(enrichment.x),y=log2(enrichment.y)))+
    geom_point()

  ggplot(ee,aes(x=(enrichment.x),y=(enrichment.y),color=-log10(combined_pvalue)))+
    geom_point()+
    scale_color_gradient2(low="black",mid = "orange",high = "red")

  cor.test((ees$enrichment.x),(ees$enrichment.y), method = "pearson")

  ggplot(dfx,aes(x=dfx$Expression,dfx$mean))+
    geom_point()

  head(ee[order(-ee$enrichment.y),])

  cor.test(log2(dfx$Expression),log2(dfx$mean))
}

identify_significantly_increased_in_remea <- function(df, fold.cutoff=0,qval.cutoff=0.05){
  ss <- sd(df$delta_rank)
  df.up <- subset(df, df$qvalue<qval.cutoff & df$delta_rank>ss)
  return(unique(df.up$perturbagens))
}

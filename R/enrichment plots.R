


plot_signature_markers <- function(marker_type,tumour_type,perturbagen,protein_data){


  signature_list <- ReMEA::list_of_signatures

  names(signature_list)

  ss <- signature_list$crispr_gerdes_haem_signatures
  # Get signatues to analyse
  signatures_to_analyse <- foreach(m=marker_type, .combine="c")%do%{
    foreach (tt = tumour_type, .combine = "c")%do%{
      return(names(signature_list)[grepl(m,names(signature_list),fixed = T) &
                                     grepl(tt,names(signature_list),fixed = T)])
    }
  }

  list.of.plots.1 <- list()
  list.of.plots.2 <- list()
  list.of.plots.3 <- list()
  list.of.plots.4 <- list()
  i <- 1
  for (signature_name in signatures_to_analyse){

    x <- ReMEA::plot_one_signature_markers(signature_name = signature_name,
                                           perturbagen = perturbagen,
                                           protein_data = protein_data)

    list.of.plots.1[[i]] <- x$plot_all
    list.of.plots.2[[i]] <- x$plot_rank
    list.of.plots.3[[i]] <- x$plot_ecdf
    list.of.plots.4[[i]] <- x$plot_density
    #x$caption
    i <- i+1
  }

  names(list.of.plots.1) <- signatures_to_analyse
  names(list.of.plots.2) <- signatures_to_analyse
  names(list.of.plots.3) <- signatures_to_analyse
  names(list.of.plots.4) <- signatures_to_analyse

  return(list(plot_all=list.of.plots.1,
              plot_rank=list.of.plots.2,
              plot_ecdf=list.of.plots.3,
              plot_density=list.of.plots.4,
              caption=x$caption))

}


plot_one_signature_markers <- function(signature_name,
                                       perturbagen,
                                       protein_data){


  colnames(protein_data)[2] <- "fold"

  protein_data <- na.omit(protein_data)
  protein_data[,2] <- scale(protein_data[,2])
  protein_data$protein_rank <- rank(protein_data[,2]) ##### Rank

  signature_list <- ReMEA::list_of_signatures


  df.sign <- signature_list[[signature_name]]
  xx <- df.sign[df.sign$perturbagen==perturbagen,]

  xx.sensitivity <- unlist(strsplit(xx$proteins.sensitivity.markers,";"))
  xx.resistance <- unlist(strsplit(xx$proteins.resistance.markers,";"))

  df.sen <- protein_data[protein_data$protein_acc %in% xx.sensitivity,]
  df.res <- protein_data[protein_data$protein_acc %in% xx.resistance,]
  df.other <- protein_data[!(protein_data$protein_acc %in% xx.resistance) &
                             !(protein_data$protein_acc %in% xx.sensitivity),]

  df.sen$signature_type <- "Sensitivity"
  df.res$signature_type <- "Resistance"
  df.other$signature_type <- "None"

  df.all <- rbind.data.frame(df.sen,df.res,df.other)

  tt <- ks.test(df.res$fold,df.sen$fold)
  ks.stat <- round(tt$statistic, digits = 2)
  ks.pvalue <- formatC(tt$p.value,format = 'e',digits = 2)

  av.rank.sen <- median(df.sen$protein_rank)
  av.rank.res <- median(df.res$protein_rank)
  delta.rank <- av.rank.sen-av.rank.res

  pp2 <- ggplot()+
    geom_point(data=df.other, aes(x=protein_rank, y=fold),color="grey",
               shape=16,alpha=0.5)+
    geom_point(data=df.sen, aes(x=protein_rank, y=fold,color="Sensitivity"),
               shape=16)+
    geom_point(data=df.res, aes(x=protein_rank, y=fold,color="Resistance"),
               shape=16)+
    scale_color_manual(values = c("red","royalblue"))+
    theme_bw()+
    labs(title = perturbagen, subtitle = signature_name, colour="signature_type",y="Expression")+
    geom_hline(yintercept = 0, linetype=2)+
    annotate(geom = 'text',y=max(protein_data$fold),x=median(protein_data$protein_rank),
             label=paste0("Median Rank Sensitivity = ", av.rank.sen, ", n = ", nrow(df.sen)))+
    annotate(geom = 'text',y=max(protein_data$fold)-1,x=median(protein_data$protein_rank),
             label=paste0("Median Rank Resistance = ", av.rank.res,", n = ", nrow(df.res))) +
    annotate(geom = 'text',y=max(protein_data$fold)-2,x=median(protein_data$protein_rank),
             label=paste0("Delta Rank (res-sen) = ", delta.rank))+
    annotate(geom='text',y=max(protein_data$fold)-3,x=median(protein_data$protein_rank),label=paste0('KS p-value = ', ks.pvalue))

  pp2

  plot.rank <- ggplot(data=rbind.data.frame(df.res,df.sen),
                      aes(x=signature_type,y=protein_rank,
                          color=signature_type))+
    # geom_boxplot()+
    geom_violin()+
  geom_point()+
    scale_color_manual(values = c("red","royalblue"))+
    theme_bw()+
    labs(title = perturbagen,
         subtitle = signature_name, fill="signature_type",y="Expression")


  plot.rank

  pp10 <- ggplot()+
    geom_point(data=df.other, aes(x=protein_rank, y=fold),color="grey",
               shape=16,alpha=0.5)+

    geom_bar(data=df.sen, aes(x=protein_rank, y=fold,fill="Sensitivity"),
             stat='identity')+
    geom_bar(data=df.res, aes(x=protein_rank, y=fold,fill="Resistance"),
             stat='identity')+
    scale_fill_manual(values = c("red","royalblue"))+
    theme_bw()+
    labs(title = perturbagen,
         subtitle = signature_name, fill="signature_type",y="Expression")+
    geom_hline(yintercept = 0, linetype=2)+
    annotate(geom = 'text',y=max(protein_data$fold),x=median(protein_data$protein_rank),
             label=paste0("Median Rank Sensitivity = ", av.rank.sen, ", n = ", nrow(df.sen)))+
    annotate(geom = 'text',y=max(protein_data$fold)-1,x=median(protein_data$protein_rank),
             label=paste0("Median Rank Resistance = ", av.rank.res,", n = ", nrow(df.res))) +
    annotate(geom = 'text',y=max(protein_data$fold)-2,x=median(protein_data$protein_rank),
             label=paste0("Delta Rank (res-sen) = ", delta.rank))+
    annotate(geom='text',y=max(protein_data$fold)-3,x=median(protein_data$protein_rank),label=paste0('KS p-value = ', ks.pvalue))

  pp10


  pp3 <- ggplot(data=rbind.data.frame(df.res,df.sen),
                aes(x=fold))+
    stat_ecdf(aes(color=signature_type), geom = 'step')+
    scale_color_manual(values = c("red","royalblue"))+
    theme_bw()+
    labs(title = perturbagen, subtitle = signature_name, y = "Fn(x)")+
    geom_hline(yintercept = 0.5, linetype=2)+
    geom_hline(yintercept = 0.25, linetype=3)+
    geom_hline(yintercept = 0.75, linetype=3)+
    annotate(geom='text',x=min(protein_data$fold)+4,y=0.95,label=paste0('KS D = ', ks.stat))+
    annotate(geom='text',x=min(protein_data$fold)+4,y=0.85,label=paste0('KS p-value = ', ks.pvalue))


  pp3



  pp4 <- ggplot(data= rbind.data.frame(df.res,df.sen),
                aes(x=signature_type,y=fold, fill=signature_type))+
    #geom_histogram(aes(y=after_stat( density)), alpha=0.5,
  #                 position="identity",bins=30)+
  geom_boxplot()+
  scale_fill_manual(values = c("red","royalblue"))+
  #scale_color_manual(values = c("grey","red","royalblue"))+
  scale_x_discrete(limits=c("Resistance","None","Sensitivity"))+
  theme_bw()+
    #ggpubr::stat_compare_means()+
    labs(title = perturbagen,
         subtitle = paste(signature_name,"\n","KS, p =", ks.pvalue))#+
    #geom_density(alpha=.2)
  pp4
  ll <- get_legend(pp4)
  px1 <- cowplot::plot_grid(pp10+theme(legend.position = "none"),
                            cowplot::plot_grid(pp3+theme(legend.position = "none"),
                            pp4+theme(legend.position = "none"),
                            ll, rel_widths = c(1,1,0.4), nrow=1),
                            nrow = 2)


  mycaption <- 'a, Proteins in resistance (red dots) and sensitivity signatures (blue dots) ranked by expression
                    b, Cummulative frequency distrubitions of proteins in resistance (red dots) and sensitivity signatures
                    c, Density distributions of proteins in resistance (red dots) and sensitivity signatures'

  pp5 <- ggplot()+annotate(geom = 'text',x=0,y=0,
                           label='a, Proteins in resistance (red dots) and sensitivity signatures (blue dots) ranked by expression
                    b, Cummulative frequency distrubitions of proteins in resistance (red dots) and sensitivity signatures
                    c, Density distributions of proteins in resistance (red dots) and sensitivity signatures',
                           hjust = 0.5)+
    theme_void()

  pp6 <- cowplot::plot_grid(px1,pp5,nrow = 2,
                            rel_heights = c(1,0.3))

  return(list(plot_all=px1,
              plot_rank =plot.rank,
              plot_ecdf=pp3,
              plot_density=pp4,
              caption=mycaption))
}




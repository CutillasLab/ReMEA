


plots_from_ReMEA <- function(df.remea_data,
                             df.scores,
                              graph_title="",
                              plot_top_n=50,
                             qvalue_cutoff=0.05){


  df.scores <- df.scores[ df.scores$qvalue<qvalue_cutoff & df.scores$av.effect>0.1,]
  df.scores <- df.scores[order(-df.scores$av.effect),]
  if (nrow(df.scores)>plot_top_n){

    df.scores <- df.scores[1:plot_top_n,]
  }
  df.scores <- df.scores[order(-df.scores$av.effect),]
  top_pers <- df.scores$perturbagen

  mm <- mean(df.remea_data$max_delta_score)
  ss <- sd(df.remea_data$max_delta_score)


  dfxx <- df.remea_data[df.remea_data$perturbagens %in% top_pers,]
  dfxx <-dfxx[order(-dfxx$max_delta_score),]
  dfxx$perturbagens <-  stringr::str_trunc(dfxx$perturbagens,30)
  pp3 <- ggplot(dfxx,aes(x=max_delta_score,y=perturbagens))+
    #geom_boxplot()+
    geom_point(aes(size=-log10(qvalue),color=-log10(qvalue)))+#, shape=signature.type))+
    #scale_shape_manual(values = c(1:20,1:20))+
    scale_color_gradient2(low="white",mid = "orange",high = "purple4",midpoint = 0.7)+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by ReMEA score",
         x="ReMEA score")+
    #theme_bw()+
    #theme(legend.position="none")+
    scale_y_discrete(limits=rev(df.scores$perturbagen))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = 0)#+
    #geom_vline(xintercept = mm+ss,linetype=2)+
    #geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
    #geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
    #annotate(geom = "text",y=2,x=mm+ss,
             #label=expression(sigma),size=10,hjust=0 ,color="grey")+
    #annotate(geom = "text",y=2,x=mm+(ss*2),
            # label=expression("2"* sigma),size=10,hjust=0,color="grey")
  pp3
  head(df.scores)
  pp3b <- ggplot(df.scores,aes(x=av.effect,y=(perturbagen)))+
    #geom_boxplot()+
    geom_point(aes(size=-log10(qvalue),color=-log10(qvalue)))+
    scale_y_discrete(limits=rev(df.scores$perturbagen))+
  #, shape=signature.type))+
    #scale_shape_manual(values = c(1:20,1:20))+
    scale_color_gradient2(low="white",mid = "orange",high = "purple4",midpoint = 0.7)+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Combined adjusted pvalue",
         color="Combined qval")#+
    #theme_bw()
    #theme(legend.position="bottom")#+
    #theme_void()
  pp3b
  pp3
  llegend <- get_legend(pp3b)
  pplot_combined <- cowplot::plot_grid(pp3+ theme_classic()+ theme(legend.position="none"),
                     pp3b+ theme_classic()+ theme(legend.position="none"),llegend, nrow = 1, rel_widths = c(1,0.5,0.3))

  return(list(plot_by_qvalue=pp3b,
              plot_by_score=pp3,
              #plot_by_deltarank=pp3,
              combined_plot=pplot_combined))


}

plot_remea_of_genes_targeted_by_drugs <- function(df.remea_data,
                                         df.scores,
                                         graph_title="",
                                         plot_top_n=50,
                                         qvalue_cutoff=0.25){


  x <- ReMEA::drug_info_from_Shelleckchem_pharmacoDB_plus_Finland$targets.genes
  genes <- unlist(strsplit(unlist(x),";"))

  df.scores <- df.scores[df.scores$perturbagen %in% genes,]
  df.remea_data <- df.remea_data[df.remea_data$perturbagens %in% genes,]



  df.scores <- df.scores[ df.scores$qvalue<qvalue_cutoff & df.scores$av.effect>0.1,]
  df.scores <- df.scores[order(-df.scores$av.effect),]
  if (nrow(df.scores)>plot_top_n){

    df.scores <- df.scores[1:plot_top_n,]
  }
  df.scores <- df.scores[order(-df.scores$av.effect),]
  top_pers <- df.scores$perturbagen

  mm <- mean(df.remea_data$max_delta_score)
  ss <- sd(df.remea_data$max_delta_score)


  dfxx <- df.remea_data[df.remea_data$perturbagens %in% top_pers,]
  dfxx <-dfxx[order(-dfxx$max_delta_score),]
  dfxx$perturbagens <-  stringr::str_trunc(dfxx$perturbagens,30)
  pp3 <- ggplot(dfxx,aes(x=max_delta_score,y=perturbagens))+
    #geom_boxplot()+
    geom_point(aes(size=-log10(qvalue),color=-log10(qvalue)))+#, shape=signature.type))+
    #scale_shape_manual(values = c(1:20,1:20))+
    scale_color_gradient2(low="white",mid = "orange",high = "purple4",midpoint = 0.7)+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by ReMEA score",
         x="ReMEA score")+
    #theme_bw()+
    #theme(legend.position="none")+
    scale_y_discrete(limits=rev(df.scores$perturbagen))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = 0)#+
  #geom_vline(xintercept = mm+ss,linetype=2)+
  #geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
  #geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
  #annotate(geom = "text",y=2,x=mm+ss,
  #label=expression(sigma),size=10,hjust=0 ,color="grey")+
  #annotate(geom = "text",y=2,x=mm+(ss*2),
  # label=expression("2"* sigma),size=10,hjust=0,color="grey")
  pp3
  head(df.scores)
  pp3b <- ggplot(df.scores,aes(x=av.effect,y=(perturbagen)))+
    #geom_boxplot()+
    geom_point(aes(size=-log10(qvalue),color=-log10(qvalue)))+
    scale_y_discrete(limits=rev(df.scores$perturbagen))+
    #, shape=signature.type))+
    #scale_shape_manual(values = c(1:20,1:20))+
    scale_color_gradient2(low="white",mid = "orange",high = "purple4",midpoint = 0.7)+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Combined adjusted pvalue",
         color="Combined qval")#+
  #theme_bw()
  #theme(legend.position="bottom")#+
  #theme_void()
  pp3b
  pp3
  llegend <- get_legend(pp3b)
  pplot_combined <- cowplot::plot_grid(pp3+ theme_classic()+ theme(legend.position="none"),
                                       pp3b+ theme_classic()+ theme(legend.position="none"),llegend, nrow = 1, rel_widths = c(1,0.5,0.3))

  return(list(plot_by_qvalue=pp3b,
              plot_by_score=pp3,
              #plot_by_deltarank=pp3,
              combined_plot=pplot_combined))

  }




plots_from_ReMEA_selected <- function(df, selected_perturbagens, graph_title=""){


  dfxx <- df[df$perturbagens %in% selected_perturbagens,]
  dfxx <-dfxx[order(-dfxx$delta_rank_median),]
  dfxx$perturbagens <-  stringr::str_trunc(dfxx$perturbagens,30)
  pp.rnai.s <- ggplot(dfxx,aes(x=-log10(qvalue),y=reorder(perturbagens, delta_rank)))+
    geom_point(aes(size=total_counts,color=-log10(qvalue)))+
    scale_color_gradient2(low="yellow",mid = "orange",high = "purple4")+
    labs(title = graph_title,
         subtitle = "p-values of response by delta rank",y="Perturbagen")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    theme_bw()

  mm <- mean(df$delta.zscore)
  ss <- sd(df$delta.zscore)

  dfxx <-dfxx[order(-dfxx$delta.zscore),]
  pp2 <- ggplot(dfxx,aes(x=delta.zscore,y=reorder(perturbagens, delta.zscore)))+
    geom_boxplot()+
    geom_point(aes(size=total_counts,color=-log10(qvalue)))+
    scale_color_gradient2(low="yellow",mid = "orange",high = "purple4")+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by delta counts")+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = mm)+
    geom_vline(xintercept = mm+ss,linetype=2)+
    geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
    geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
    annotate(geom = "text",y=2,x=mm+ss,
             label=expression(sigma),size=10,hjust=0 ,color="grey")+
    annotate(geom = "text",y=2,x=mm+(ss*2),
             label=expression("2"* sigma),size=10,hjust=0,color="grey")

  pp2

  mm <- mean(df$delta_rank_median)
  ss <- sd(df$delta_rank_median)


  dfxx <-dfxx[order(-dfxx$delta_rank_median),]
  pp3.med <- ggplot(dfxx,aes(x=delta_rank_median,y=reorder(perturbagens, delta_rank_median)))+
    geom_boxplot()+
    geom_point(aes(size=total_counts,color=-log10(qvalue)))+
    scale_color_gradient2(low="yellow",mid = "orange",high = "purple4")+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by delta rank")+theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = mm)+
    geom_vline(xintercept = mm+ss,linetype=2)+
    geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
    geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
    annotate(geom = "text",y=2,x=mm+ss,
             label=expression(sigma),size=10,hjust=0 ,color="grey")+
    annotate(geom = "text",y=2,x=mm+(ss*2),
             label=expression("2"* sigma),size=10,hjust=0,color="grey")
  pp3.med
  mm <- mean(df$delta_rank_mean)
  ss <- sd(df$delta_rank_mean)

  dfxx <-dfxx[order(-dfxx$delta_rank_mean),]
  pp3.mean <- ggplot(dfxx,aes(x=delta_rank_mean,y=reorder(perturbagens, delta_rank_mean)))+
    geom_boxplot()+

    geom_point(aes(size=total_counts,color=-log10(qvalue)))+
    scale_color_gradient2(low="yellow",mid = "orange",high = "purple4")+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by delta rank")+theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = mm)+
    geom_vline(xintercept = mm+ss,linetype=2)+
    geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
    geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
    annotate(geom = "text",y=2,x=mm+ss,
             label=expression(sigma),size=10,hjust=0 ,color="grey")+
    annotate(geom = "text",y=2,x=mm+(ss*2),
             label=expression("2"* sigma),size=10,hjust=0,color="grey")
  pp3.mean

  mm <- mean(df$delta_rank_geommean)
  ss <- sd(df$delta_rank_geommean)

  dfxx <-dfxx[order(-dfxx$delta_rank_geommean),]
  pp3.geommean <- ggplot(dfxx,aes(x=delta_rank_geommean,y=reorder(perturbagens, delta_rank_geommean)))+
    geom_boxplot()+
    geom_point(aes(size=total_counts,color=-log10(qvalue)))+
    scale_color_gradient2(low="yellow",mid = "orange",high = "purple4")+
    labs(title = graph_title , y="perturbagens",
         subtitle = "Predicted perturbagens response, ranked by delta rank")+theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))+
    geom_vline(xintercept = mm)+
    geom_vline(xintercept = mm+ss,linetype=2)+
    geom_vline(xintercept = mm+ss,linetype=2,color="grey")+
    geom_vline(xintercept = mm+(2*ss),linetype=2,color="grey")+
    annotate(geom = "text",y=2,x=mm+ss,
             label=expression(sigma),size=10,hjust=0 ,color="grey")+
    annotate(geom = "text",y=2,x=mm+(ss*2),
             label=expression("2"* sigma),size=10,hjust=0,color="grey")
pp3.geommean
  ppx.s <- cowplot::plot_grid(pp2, pp3.med, pp3.mean,pp3.geommean,
                              rel_widths = c(1,1))
  ppx.s
  return (list(plot_combined=ppx.s,
             plot_by_deltarank=pp2.rnai.s,
             plot_by_qvalue=pp.rnai.s))
}

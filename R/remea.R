




remea <- function(protein_data, tumour_type=c("haem","solid"), expt_title=""){


  SigEnrichDrug <-  ReMEA::response_marker_enrichment_analysis_v2(protein_data = na.omit(protein_data),
                                                                  marker_type = "drug",
                                                                  tumour_type = tumour_type)


  drug_score_data <-  ReMEA::combined_remea_score(df.remea_results=SigEnrichDrug$enrichment_data)

  plot_remea_drug <- drug_score_data$remea_plot+ggtitle("Drug ReMEA",subtitle = expt_title)





  SigEnrichCrisp <-  ReMEA::response_marker_enrichment_analysis_v2(protein_data = na.omit(protein_data),
                                                                   marker_type = c("rnai","crispr"),
                                                                   tumour_type = tumour_type)

  gene_score_data <-  ReMEA::combined_remea_score(df.remea_results=SigEnrichCrisp$enrichment_data)
  plot_remea_gene <- gene_score_data$remea_plot+ggtitle("Gene silencing ReMEA",subtitle = expt_title)
  plot_remea_gene <- drug_score_data$remea_plot+ggtitle("Drug ReMEA",subtitle = expt_title)



  x <- gene_score_data$combined_score
  x[x$perturbagen=="BCL2",]

  ###########################

  drug_target_genes <- unique(unlist(strsplit( SigEnrichDrug$enrichment_data$targets.genes,";")))
  x2 <- x[x$perturbagen %in% drug_target_genes,]
  x2 <- subset(x2,x2$qvalue<0.05)
  nrow(x2)
  x2 <- x2[order(-x2$alpha),]

  if (nrow(x2)>50){
    x2 <- x2[1:50,]
  }

  pplot <- ggplot(x2,aes(x=av.effect,y=reorder(perturbagen,av.effect),color=-log10(qvalue)))+
    geom_point(size=5)+
    scale_color_gradient2(low="grey",mid = "orange",high = "red", midpoint = median(-log10(x2$qvalue),na.rm = T))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    ggtitle("ReMEA of drug target genes ",subtitle = expt_title)

  pplot
  #############################


  x <- gene_score_data$combined_score
  x4 <- subset(x,x$qvalue<0.05 & x$av.effect>0.02)
  nrow(x4)
  x4 <- x4[order(-x4$alpha),]
  if (nrow(x4)>50){
    x4 <- x4[1:50,]
  }
  pplot_genes <- ggplot(x4,aes(x=av.effect,y=reorder(perturbagen,av.effect),color=-log10(qvalue)))+
    geom_point(size=5)+
    scale_color_gradient2(low="grey",mid = "orange",high = "red", midpoint = median(-log10(x2$qvalue),na.rm = T))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    ggtitle("ReMEA of drug target genes ",subtitle = expt_title)

  pplot_genes
  #######################

  x3 <- drug_score_data$combined_score
  x3 <- subset(x3,x3$qvalue<0.05 & x3$av.effect>0.2)
  x3 <- x3[order(-x3$alpha),]
  x3$perturbagen <-  stringr::str_trunc(x3$perturbagen,30)
  if (nrow(x3)>50){
    x3 <- x3[1:50,]
  }
  pplot_drug <- ggplot(x3,aes(x=av.effect,y=reorder(perturbagen,av.effect),color=-log10(qvalue)))+
    geom_point(size=5)+
    scale_color_gradient2(low="grey",mid = "orange",high = "red",
                          midpoint = median(-log10(x2$qvalue),na.rm = T))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    ggtitle("ReMEA of drug target genes ",subtitle = expt_title)
  pplot_drug
  x3s <- x3[1:30,]
  pplot_drug_all <- ggplot(x3,aes(x=av.effect,y=-log10(qvalue)))+
    geom_point(aes(size=alpha,color=alpha))+
    geom_text_repel(data = x3s,aes(x=av.effect,y=-log10(qvalue),label=perturbagen))+
    scale_color_gradient2(low="grey",mid = "orange",high = "red", midpoint = median(-log10(x2$qvalue),na.rm = T))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    ggtitle("ReMEA of drug target genes ",subtitle = expt_title)

  pplot_drug_all


  return(list(remea_scores_drug_combined= drug_score_data$combined_score,
              remea_scores_gene_combined= gene_score_data$combined_score,
              remea_scores_drug_individual_dbs = SigEnrichDrug$enrichment_data,
              remea_scores_gene_individual_dbs = SigEnrichCrisp$enrichment_data,
              plot_drug_scores = pplot_drug,
              plot_gene_scores = pplot_genes,
              plot_drug_target_genes_scores = pplot

  ))


  }



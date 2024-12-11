#' Plot combined ReMEA score
#'
#' @param combined_remea `data.table` of combined ReMEA scores. Output of `combine_remea_scores` function.
#' @param drug_target_genes `logical` Whether to plot only drug target genes.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_ReMEA_scores <- function(combined_remea,
                              drug_target_genes = FALSE){
  if (!data.table::is.data.table(combined_remea)) {
    message("Note, combined_remea was not passed as data.table, converting now...")
    combined_remea <- data.table::as.data.table(combined_remea)
  }
  if (drug_target_genes == TRUE){
    subset_dt <- combined_remea[order(-av.effect)][perturbagen %in% ReMEA::drug_target_genes][1:50]
  } else {
    subset_dt <- combined_remea[order(-av.effect)][1:50]
  }
  # Replace 0 p values
  cols_to_replace <- c("ks_qvalue", "bws_qvalue")
  subset_dt[, (cols_to_replace) := lapply(.SD, function(x) ifelse(x == 0, 1e-300, x)), .SDcols = cols_to_replace]
  subset_dt[, alpha := av.effect * -log10(ks_qvalue)]
  ggplot(subset_dt,aes(x=av.effect, y=reorder(perturbagen, av.effect),color=-log10(ks_qvalue)))+
    geom_point(aes(size=alpha))+
    scale_color_gradient2(low="grey",mid = "orange",high = "red", midpoint = median(-log10(subset_dt$ks_qvalue),na.rm = T))+
    theme_bw()+
    theme(axis.title.y = element_blank(), legend.direction = "vertical")+
    xlab("Combined ReMEA score")+
    geom_vline(xintercept = 0)+
    geom_vline(xintercept = 1, linetype = "dashed") -> remea_plot
  return(remea_plot)
}

#' Plot individual database max delta score
#'
#' @param remea_db_scores `data.table` of individual db scores. Output of `response_marker_enrichment_analysis` function.
#' @param drug_target_gens `logical` Whether to plot only drug target genes.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_individual_db_ReMEA_scores <- function(remea_db_scores,
                                            drug_target_genes = FALSE){
  # Check if object is data.table.
  if (!data.table::is.data.table(remea_db_scores)) {
    print("Note, remea_results was not passed as data.table, converting now...")
    remea_db_scores <- data.table::as.data.table(remea_db_scores)
  }
  if (drug_target_genes == TRUE){
    top_50 <- remea_db_scores[perturbagen %in% ReMEA::drug_target_genes,
                         .(mean_mds = mean(max_delta_score)),
                         by = perturbagen][order(-mean_mds)][1:50][,perturbagen]
  } else {
    top_50 <- remea_db_scores[, .(mean_mds = mean(max_delta_score)), by = perturbagen][order(-mean_mds)][1:50][,perturbagen]
  }
  # P value color
  min_p <- -log10(min(remea_db_scores[perturbagen %in% top_50, ks_qvalue]))
  max_p <- -log10(max(remea_db_scores[perturbagen %in% top_50, ks_qvalue]))
  if (is.infinite(min_p )){min_p = 0}
  points_p <- seq(max_p, min_p, length.out = 4)

  ggplot(remea_db_scores[perturbagen %in% top_50],
         aes(x = max_delta_score, y = reorder(perturbagen, max_delta_score)))+
    geom_point(aes(size = -log10(ks_qvalue), color = -log10(ks_qvalue)))+
    scale_color_gradient2(low="white", mid = "orange", high = "purple4", midpoint = 0.7, breaks=points_p)+
    scale_size_continuous(breaks=points_p)+
    theme_bw()+
    guides(color= guide_legend(), size=guide_legend())+
    theme(axis.title.y = element_blank(),
          legend.direction = "vertical")+
    xlab("Single DB max delta scores")+
    geom_vline(xintercept = 0)+
    geom_vline(xintercept = 1,
               linetype = "dashed") -> remea_db_plot
  return(remea_db_plot)
}

#' Plot ReMEA score volcano
#'
#' @param combined_remea `data.table` of combined ReMEA scores. Output of `combine_remea_scores` function.
#' @param qval_cut `float` labelling cut off for volcano plot.
#' @param drug_target_genes `logical` whether to just plot all genes or just drug target genes.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_ReMEA_volcano <- function(combined_remea,
                               qval_cut = 0.05,
                               drug_target_genes = FALSE){
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("This function requires the ggrepel library")
  }
  # Check if object is data.table.
  if (!data.table::is.data.table(combined_remea)) {
    message("Note, 'combined_remea' was not passed as data.table, converting now...")
    combined_remea <- data.table::as.data.table(combined_remea)
  }
  if (drug_target_genes == TRUE){
    combined_remea <- combined_remea[perturbagen %in% ReMEA::drug_target_genes]
  }
  ggplot(combined_remea,aes(x=av.effect,y=-log10(ks_qvalue)))+
    geom_point(aes(size=abs(alpha),color=alpha))+
    ggrepel::geom_text_repel(data = combined_remea[ks_qvalue < qval_cut],aes(x=av.effect,y=-log10(ks_qvalue),label=perturbagen))+
    scale_color_gradient2(low="blue",mid = "orange",high = "red", midpoint = median(-log10(combined_remea$ks_qvalue),na.rm = T))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    geom_vline(xintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = -1, linetype = "dashed") -> volcano_plot
  return(volcano_plot)
}

#' Plot CRISPR versus RNAi ReMEA scores
#'
#' @param remea_scores `list` output from `complete_ReMEA_analysis` function.
#' @param combined_RNAi `data.table` of combined RNAi ReMEA scores.
#' @param combined_CRISPR `data.table` of combined CRISPR ReMEA scores.
#' @param drug_target_genes `logical` whether to filter plot to just include genes of drug targets.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_crispr_vs_rnai <- function(remea_scores = NULL,
                                combined_RNAi = NULL,
                                combined_CRISPR = NULL,
                                drug_target_genes = FALSE){
  if (is.null(remea_scores) & is.null(combined_RNAi) & is.null(combined_CRISPR)){
    stop("Please pass either: List of complete ReMEA analysis, or combined scores for RNAi and CRISPR perturbations.")
  } else if (!is.null(remea_scores)){
    combined_RNAi <- remea_scores$remea_scores_rnai_combined
    combined_CRISPR <- remea_scores$remea_scores_crispr_combined
  } else if (is.null(combined_RNAi) | is.null(combined_CRISPR)){
    stop("Please pass the required data.")
  }
  combined_RNAi <- combined_RNAi[, .(perturbagen, av.effect, ks_qvalue, alpha)]
  colnames(combined_RNAi)[-1] <- paste0("RNAi_", colnames(combined_RNAi)[-1])
  combined_CRISPR <- combined_CRISPR[, .(perturbagen, av.effect, ks_qvalue, alpha)]
  colnames(combined_CRISPR)[-1] <- paste0("CRISPR_", colnames(combined_CRISPR)[-1])
  rnai_crispr <- data.table::merge.data.table(combined_RNAi, combined_CRISPR)
  if (drug_target_genes == TRUE){
    rnai_crispr <- rnai_crispr[perturbagen %in% ReMEA::drug_target_genes]
  }
  cols_to_replace <- c("RNAi_ks_qvalue", "CRISPR_ks_qvalue")
  rnai_crispr[, (cols_to_replace) := lapply(.SD, function(x) ifelse(x == 0, 1e-300, x)), .SDcols = cols_to_replace]
  rnai_crispr[, `:=`(RNAi_alpha = RNAi_av.effect * -log10(RNAi_ks_qvalue),
                     CRISPR_alpha = CRISPR_av.effect * -log10(CRISPR_ks_qvalue))]
  rnai_crispr[, avg.alpha := rowMeans(.SD), .SDcols = c("RNAi_alpha", "CRISPR_alpha")]
  limit_max <- max(max(rnai_crispr[, CRISPR_av.effect]),
                   max(rnai_crispr[, RNAi_av.effect]))
  limit_min <- min(min(rnai_crispr[, CRISPR_av.effect]),
                   min(rnai_crispr[, RNAi_av.effect]))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggplot(data = rnai_crispr, aes(x = RNAi_av.effect, CRISPR_av.effect, label = perturbagen))+
      geom_rect(aes(xmin=1,xmax=Inf,ymin=1,ymax=Inf),alpha=0.05,fill="#E6E6E3")+
      geom_point(aes(color = avg.alpha))+
      theme_bw()+
      xlim(c(limit_min - 0.5, limit_max + 0.5))+
      ylim(c(limit_min - 0.5, limit_max + 0.5))+
      scale_color_gradient2(low="blue",mid = "orange",high = "red", midpoint = median(rnai_crispr$avg.alpha,na.rm = T))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_hline(yintercept = c(-1,1), color = "#6F7378")+
      geom_vline(xintercept = c(-1,1), color = "#6F7378")+
      ggrepel::geom_text_repel(data = rnai_crispr[(CRISPR_av.effect > 1 & RNAi_av.effect > 1) | (CRISPR_av.effect < -1 & RNAi_av.effect < -1)],
                               max.overlaps = 50, force =3,
                               alpha = 0.6)+
      theme(axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 14))-> remea_plot
  } else {
    stop("The 'ggrepel' package is required for this plot.")
  }
  return(remea_plot)
}

#' Plot drug ReMEA scores versus gene target perturbation ReMEA scores
#'
#' @param remea_scores `list` output from `complete_ReMEA_analysis` function.
#' @param combined_drug `data.table` of combined RNAi ReMEA scores.
#' @param combined_gene `data.table` of combined CRISPR ReMEA scores.
#' @param gene_perturbation_type `string` indicating whether "RNAi", "CRISPR" or "RNAi_CRISPR" scores should be extracted from `remea_scores` list.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_drug_vs_gene <- function(remea_scores = NULL,
                              combined_drug = NULL,
                              combined_gene = NULL,
                              gene_perturbation_type = NULL){
  if (is.null(remea_scores) & is.null(combined_drug) & is.null(combined_gene)){
    print("Please pass either: List of complete ReMEA analysis, or combined scores for RNAi and CRISPR perturbations.")
  } else if (!is.null(remea_scores)){
    if (is.null(gene_perturbation_type)){
      print("If passing complete_ReMEA_analysis list, please choose genetic perturbation type.")
    } else if (gene_perturbation_type == "RNAi"){
      combined_gene <- remea_scores$remea_scores_rnai_combined
    } else if (gene_perturbation_type == "CRISPR"){
      combined_gene <- remea_scores$remea_scores_crispr_combined
    } else if (gene_perturbation_type == "RNAi_CRISPR") {
      combined_gene <- remea_scores$remea_scores_gene_combined
    } else {
      print("Please supply correct input.")
    }
    combined_drug <- remea_scores$remea_scores_drug_combined
  } else if (is.null(combined_drug) | is.null(combined_gene)){
    print("Please pass the required data.")
  }

  target_drug_map <- as.data.table(ReMEA::drug_target_map)
  target_drug_map[, drug_target := paste0(drug_name, "__", gene_target)]
  drug_temp<-data.table::merge.data.table(target_drug_map, combined_drug[, .(perturbagen, av.effect, ks_qvalue)], by.x = "drug_name", by.y = "perturbagen")
  drug_temp <- drug_temp[!is.na(gene_target)]
  gene_temp<-data.table::merge.data.table(target_drug_map, combined_gene[, .(perturbagen, av.effect, ks_qvalue)], by.x = "gene_target", by.y = "perturbagen")
  setnames(drug_temp, c("av.effect", "ks_qvalue"), c("DRUG_av.effect", "DRUG_ks_qvalue"))
  setnames(gene_temp, c("av.effect", "ks_qvalue"), c("GENE_av.effect", "GENE_ks_qvalue"))

  drug_target_remeas <- data.table::merge.data.table(x = drug_temp[, .(drug_target, drug_name, DRUG_av.effect, DRUG_ks_qvalue, pathway)],
                                                     y = gene_temp[, .(drug_target , gene_target, GENE_av.effect, GENE_ks_qvalue)],
                                                     by = "drug_target")
  cols_to_replace <- c("DRUG_ks_qvalue", "GENE_ks_qvalue")
  drug_target_remeas[, (cols_to_replace) := lapply(.SD, function(x) ifelse(x == 0, 1e-300, x)), .SDcols = cols_to_replace]
  drug_target_remeas[, `:=`(DRUG_alpha = DRUG_av.effect * -log10(DRUG_ks_qvalue),
                            GENE_alpha = GENE_av.effect * -log10(GENE_ks_qvalue))]
  drug_target_remeas[, avg.alpha := rowMeans(.SD), .SDcols = c("DRUG_alpha", "GENE_alpha")]

  limit_max <- max(max(drug_target_remeas[, GENE_av.effect]),
                   max(drug_target_remeas[, DRUG_av.effect]))
  limit_min <- min(min(drug_target_remeas[, GENE_av.effect]),
                   min(drug_target_remeas[, DRUG_av.effect]))

  drug_target_remeas[, `Target pathway` := ifelse(GENE_av.effect > 1 & DRUG_av.effect > 1, pathway, NA)]
  drug_target_remeas <- drug_target_remeas[!duplicated(drug_target_remeas[,.(drug_target,drug_name)])]
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggplot()+
      geom_point(data = drug_target_remeas[!is.na(`Target pathway`)], aes(x = DRUG_av.effect, GENE_av.effect, size=2, color = avg.alpha, shape=`Target pathway`))+
      geom_point(data = drug_target_remeas[is.na(`Target pathway`)], aes(x = DRUG_av.effect, GENE_av.effect,color = avg.alpha))+
      theme_bw()+
      geom_rect(aes(xmin=1,xmax=Inf,ymin=1,ymax=Inf),alpha=0.5,fill="#E6E6E3")+
      scale_color_gradient2(low="blue",mid = "orange",high = "red", midpoint = median(drug_target_remeas$avg.alpha))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_hline(yintercept = c(-1,1), color = "#6F7378")+
      geom_vline(xintercept = c(-1,1), color = "#6F7378")+
      ggrepel::geom_text_repel(data = drug_target_remeas[(DRUG_av.effect > 1 & GENE_av.effect > 1) | (DRUG_av.effect < -1 & GENE_av.effect < -1)],
                               aes(x = DRUG_av.effect, y = GENE_av.effect, label = drug_target),
                               max.overlaps = 50,
                               force =3,
                               box.padding = 1,
                               alpha = 0.6)+
      guides(size="none")+
      theme(axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 14))-> remea_plot
  } else {
    stop("The 'ggrepel' package is required for this plot.")
  }
  return(remea_plot)
}

#' ReMEA plot panel
#'
#' @param individual_db_scores `data.table` of individual db scores. Output of `response_marker_enrichment_analysis` function.
#' @param combined_remea `data.table` of combined ReMEA scores. Output of `combine_remea_scores` function.
#' @param drug_target_genes `data.table` of combined ReMEA scores. Output of `combine_remea_scores` function.
#'
#' @return
#' @export
remea_plots <- function(individual_db_scores,
                        combined_remea,
                        drug_target_genes = FALSE){

  if (!requireNamespace("cowplot", quietly = TRUE)){
    stop("This function requires the 'cowplot' library")
  }
  if (drug_target_genes == TRUE){
    combined_remea <- combined_remea[perturbagen %in% ReMEA::drug_target_genes]
    individual_db_scores <- individual_db_scores[perturbagen %in% ReMEA::drug_target_genes]
  }
  av.effect.plot <- ReMEA::plot_ReMEA_scores(combined_remea)
  individual.score.plot <- ReMEA::plot_individual_db_ReMEA_scores(individual_db_scores)
  volcano.plot <- ReMEA::plot_ReMEA_volcano(combined_remea)

  combined_plot <- cowplot::plot_grid(av.effect.plot,
                                      individual.score.plot,
                                      volcano.plot,
                                      ncol = 3,
                                      rel_widths = c(1.2, 1.2, 2))
}

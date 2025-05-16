#' Combine ReMEA score
#'
#' Function to average ReMEA scores for perturbagens across different datasets.
#'
#' @param remea_results Output from response_marker_enrichment_analysis_v2 function.
#'
#' @return `data.table` of combined ReMEA scores across signature DBs.
#'
#' @export
#'
#' @examples
#' combine_remea_score(remea_results = remea_results_dt)
combine_remea_scores <- function(remea_results){
  # Check if object is data.table.
  if (!data.table::is.data.table(remea_results)) {
    print("Note, remea_results was not passed as data.table, converting now...")
    remea_results <- data.table::as.data.table(remea_results)  # Convert to data.table if it's not
  }

  .scale_get_mean <- function(dtr, by="delta_rank_median"){
    xx <- data.table::dcast(perturbagen ~ signature.type, value.var = by, data = dtr)
    xx[, mean := scale(rowMeans(.SD, na.rm = TRUE)), .SDcols = is.numeric]
    return(xx)
  }

  dtx.med <- .scale_get_mean(dtr = remea_results, "delta_rank_median")
  dtx.mean <- .scale_get_mean(dtr = remea_results, "delta_rank_mean")
  dtx.geommean <- .scale_get_mean(dtr = remea_results, "delta_rank_geomean")
  dtx.zscore <- .scale_get_mean(dtr = remea_results, "delta.zscore")

  dt.combi.scores <- data.table::data.table()
  dt.combi.scores[, `:=` (perturbagen = dtx.zscore$perturbagen,
                          av.effect = rowMeans(cbind(dtx.med$mean,
                                                     dtx.mean$mean,
                                                     dtx.geommean$mean,
                                                     dtx.zscore$mean)))]

  .combine_pvalues <- function(dtr){
    xx.ks <- data.table::dcast(perturbagen ~ signature.type, value.var = "pvalue.ks",data=dtr)
    xx.ks <- xx.ks[, combined_ks_pvalue := ReMEA::pcom_by_stouffer(.SD),
                   .SDcols = is.numeric,
                   by = 1:nrow(xx.ks)][,.(perturbagen, combined_ks_pvalue)]
    xx.bws <- data.table::dcast(perturbagen ~ signature.type, value.var = "pvalue.bws",data=dtr)
    xx.bws <- xx.bws[, combined_bws_pvalue := ReMEA::pcom_by_stouffer(.SD),
                     .SDcols = is.numeric,
                     by = 1:nrow(xx.bws)][,.(perturbagen, combined_bws_pvalue)]
    return(data.table::merge.data.table(x = xx.ks,
                                        y = xx.bws,
                                        by = "perturbagen"))
  }

  dt.combi.pvalues <- .combine_pvalues(remea_results)

  dt.combi <- data.table::merge.data.table(dt.combi.scores,
                                           dt.combi.pvalues,
                                           by = "perturbagen")
  dt.combi[,
           c("alpha",
             "ks_qvalue",
             "ks_padj_bonferroni",
             "bws_qvalue",
             "bws_padj_bonferroni") := .(
               -log10(combined_ks_pvalue) * av.effect,
               p.adjust(combined_ks_pvalue, method = "fdr"),
               p.adjust(combined_ks_pvalue, method = "bonferroni"),
               p.adjust(combined_bws_pvalue, method = "fdr"),
               p.adjust(combined_bws_pvalue, method = "bonferroni")
             )]

  return(dt.combi)
}


#' pcom by stouffer
#'
#' @param p pvalue
#'
#' @return combined p values.
#' @export
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


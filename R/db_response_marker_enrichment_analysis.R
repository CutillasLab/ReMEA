#' @title Response marker enrichment analysis
#'
#' @description Performs signture delta ranks and delta zscores iteratively across
#' selected signature DBs. Determined by pertubation type, tumour type and siganture
#' set.
#' @param protein_data `data.frame` Proteomics differential testing results. Protein name
#' must be included as protein_acc. Fold change data must be in column 2.
#' @param marker_type `string` Perturbagen type for analysis. These are CRISPR, RNAi or DRUG.
#' Multiple markers may be passed. In case of selecting both RNAi and CRISPR, these
#' may be combined using the `combine_remea_score` function
#' @param tumour_type `string` Whether to use markers from correlation analysis of solid tumours
#'  or non-solid tumours only. If multiple passed, these will be analysed in looped format.
#' @param signature_version `string` Signature version to use.
#' @param dataset_selection `logical` Whether to select specific signature DBs.
#' @param proteomics_datasets `string` Ids for proteomics dataset based signature
#' DB section.
#' @param perturbation_datasets `string` Ids for perturbation dataset based signature
#' DB section.
#' @return `data.table` ReMEA scores and p values for selected signature DBs.
#' @importFrom foreach %do%
#' @import doFuture
#' @export
#'
#' @examples
#' response_marker_enrichment_analysis(protein_data,
#' marker_type=c("DRUG", "RNAi", "CRISPR"),
#' tumour_type=c("haem","solid"))
response_marker_enrichment_analysis <- function(protein_data,
                                                marker_type=c("DRUG", "RNAi", "CRISPR"),
                                                tumour_type=c("haem","solid"),
                                                signature_version = "CellLines",
                                                dataset_selection = FALSE,
                                                proteomics_datasets = NULL,
                                                perturbation_datasets = NULL){
  # Check if object is data.table.
  if (data.table::is.data.table(protein_data)) {
    protein_data <- as.data.frame(protein_data)
  }
  # Get signatures to analyse
  signature_list <- ReMEA::select_signature(signature_version = signature_version)
  signatures_to_analyse <- foreach::foreach(m=marker_type, .combine="c")%do%{
    foreach::foreach(tt = tumour_type, .combine = "c")%do%{
      return(names(signature_list)[grepl(m, names(signature_list), fixed = FALSE, ignore.case = TRUE) &
                                     grepl(tt, names(signature_list), fixed = FALSE, ignore.case = TRUE)])
    }
  }
  if (length(signatures_to_analyse) == 0){
    stop("No signature DB with those arguments has been found.
         Ensure signature id options refer to atleast one of the listed signature lists.")
  }
  # If selected, only use signatures from specific datasets
  if (dataset_selection == TRUE){
    if (!is.null(proteomics_datasets) & !is.null(perturbation_datasets)){
      sigs_prot <- grepl(paste0(proteomics_datasets, collapse = "|"), signatures_to_analyse, ignore.case = TRUE)
      sigs_pert <- grepl(paste0(perturbation_datasets, collapse = "|"), signatures_to_analyse, ignore.case = TRUE)
      signatures_to_analyse <- signatures_to_analyse[sigs_prot & sigs_pert]
    } else if (!is.null(proteomics_datasets) & is.null(perturbation_datasets)){
      sigs_prot <- grepl(paste0(proteomics_datasets,collapse = "|"), signatures_to_analyse, ignore.case = TRUE)
      signatures_to_analyse <- signatures_to_analyse[sigs_prot]
    } else if (is.null(proteomics_datasets) & !is.null(perturbation_datasets)){
      sigs_pert <- grepl(paste0(perturbation_datasets,collapse = "|"), signatures_to_analyse, ignore.case = TRUE)
      signatures_to_analyse <- signatures_to_analyse[sigs_pert]
    }
  }
  # Perform enrichment
  future::plan(future::multisession)
  result_list <- foreach::foreach(s = signatures_to_analyse,
                                  .options.future = list(seed = TRUE))%dofuture%{
    sl <- signature_list[[s]]

    if (nrow(sl)>1){
      return(ReMEA::signature_enrichment(protein_data = protein_data,
                                         signatures = sl,
                                         signature.type = s))
    }
  }
  enrichment_dt <- data.table::rbindlist(result_list)
  enrichment_dt[
    , c("delta.zscore",
        "total_counts",
        "delta_rank_mean",
        "delta_rank_median",
        "delta_rank_geomean") := .(
          zscore.sensitivity - zscore.resistance,
          n.resistance.markers + n.sensitivity.markers,
          log2(av.rank.sensitivity) - log2(av.rank.resistance),
          log2(med.rank.sensitivity) - log2(med.rank.resistance),
          log2(geommean.rank.sensitivity) - log2(geommean.rank.resistance)
        )]
  enrichment_dt[, max_delta_score := apply(.SD[, c("delta.zscore", "delta_rank_geomean",
                                                   "delta_rank_mean", "delta_rank_median")],
                                           1, function(x) max(abs(x), na.rm = T))]
  enrichment_dt[, c("ks_qvalue", "ks_p_padj_bonf", "bws_qvalue", "bws_padj_bonf") := .(
    p.adjust(pvalue.ks, method = "BH"),
    p.adjust(pvalue.ks, method = "bonferroni"),
    p.adjust(pvalue.bws, method = "BH"),
    p.adjust(pvalue.bws, method = "bonferroni")
  ), by = signature.type]

  return(enrichment_dt)
}


#' Signature enrichment
#'
#' Main marker enrichment function for ReMEA pipeline.
#'
#' @param protein_data `data.frame` Proteomics differential testing results. Protein name
#' must be included as protein_acc. Fold change data must be in column 2.
#' @param signatures `data.table` Signature DB with response protein signatures.
#' This can be found in the ReMEA NAMESPACE. e.g. ReMEA::list_of_signature_lists
#' @param signature.type `string` Naming variable for signatures. Main purpose if for
#' response_marker_enrichment_analysis, where this indicates signature DB ids.
#'
#' @return `data.table` with enrichment results
#' @export
#'
#' @examples
#' signature_enrichemnt(protein_data = protein_data, signatures  = signature_dt, signature.type = "signature_id")
signature_enrichment <- function(protein_data,
                                 signatures,
                                 signature.type){
  # Check if object is data.table.
  if (data.table::is.data.table(protein_data)) {
    protein_data <- as.data.frame(protein_data)
  }
  protein_data <- na.omit(protein_data)
  protein_data[,2] <- scale(protein_data[,2]) # How will this effect results?
  protein_data$protein_rank <- rank(protein_data[,2])
  mean.population <- mean(protein_data[,2])
  sd.population <- sd(protein_data[,2])

  dt_enrichment <- signatures[n.resistance.markers > 2 & n.sensitivity.markers > 2, {
    resistance.signature <- unlist(strsplit(proteins.resistance.markers,';'))
    sensitivity.signature <- unlist(strsplit(proteins.sensitivity.markers,';'))
    df.prot.res <- protein_data[protein_data$protein_acc %in% resistance.signature, ]
    df.prot.sen <- protein_data[protein_data$protein_acc %in% sensitivity.signature, ]

    if (nrow(df.prot.res)>2 & nrow(df.prot.sen)>2){
      resistance.markers <- paste(unique(df.prot.res$protein_acc),collapse = ";")
      sensitivity.markers <- paste(unique(df.prot.sen$protein_acc),collapse = ";")
      ks.p <- ks.test(df.prot.res[,2],df.prot.sen[,2]) # This needs investigation
      pvalue.ks <- ks.p$p.value
      pvalue.bws <- tryCatch({
        BWStest::bws_test(df.prot.res[,2], df.prot.sen[,2],
                          method = 'BWS',
                          alternative = 'two.sided')$p.value
      }, error = function(e) {
        1
      })
      zscore.resistance <- (mean(df.prot.res[,2])-mean.population)/sd.population
      zscore.sensitivity <- (mean(df.prot.sen[,2])-mean.population)/sd.population

      perturbagen.signature <- signature
      av.rank.resistance <-  mean(df.prot.res$protein_rank)
      av.rank.sensitivity <- mean(df.prot.sen$protein_rank)

      geommean.rank.resistance <- ReMEA::gm_mean(df.prot.res$protein_rank)
      geommean.rank.sensitivity <- ReMEA::gm_mean(df.prot.sen$protein_rank)

      med.rank.resistance <- median(df.prot.res$protein_rank)
      med.rank.sensitivity <- median(df.prot.sen$protein_rank)

      n.resistance.markers <- nrow(df.prot.res)
      n.sensitivity.markers <- nrow(df.prot.sen)
      .(
        zscore.resistance,
        zscore.sensitivity,
        pvalue.ks,
        pvalue.bws,
        resistance.markers,
        sensitivity.markers,
        av.rank.resistance,
        av.rank.sensitivity,
        geommean.rank.resistance,
        geommean.rank.sensitivity,
        med.rank.resistance,
        med.rank.sensitivity,
        n.resistance.markers,
        n.sensitivity.markers,
        signature.type
      )
      }},
    by = perturbagen
  ]

  return(dt_enrichment)
}

#' Get averaged ReMEA scores
#'
#' @param protein_data `data.frame` Proteomics differential testing results. Protein name
#' must be included as protein_acc. Fold change data must be in column 2.
#' @param tumour_type `string` Whether to use markers from correlation analysis of solid tumours
#' or non-solid tumours only. If multiple passed, these will be analysed in looped format.
#' @param marker_type `string` Perturbagen type for analysis. These are CRISPR, RNAi or DRUG.
#' Multiple markers may be passed. In case of selecting both RNAi and CRISPR, these
#' may be combined using the `combine_remea_score` function.
#' @param signature_version `string` Signature version to use.
#' @param return_individual_scores `logical` whether to return results for individual db scores.
#'
#' @return `data.table` of combined ReMEA scores if return_individual_scores == FALSE,
#' else returns list of both combined and individual scores.
#' @export
#'
#' @examples
get_ReMEA_scores <- function(protein_data,
                             tumour_type,
                             marker_type = "DRUG",
                             signature_version = "CellLines",
                             return_individual_scores = FALSE){
  # Check if object is data.table.
  if (data.table::is.data.table(protein_data)) {
    protein_data <- as.data.frame(protein_data)
  }
  remea_indv_dbs <- ReMEA::response_marker_enrichment_analysis(protein_data = protein_data,
                                                               marker_type = marker_type,
                                                               tumour_type = tumour_type,
                                                               signature_version = signature_version)

  combined_remea_scores <- ReMEA::combine_remea_scores(remea_results = remea_indv_dbs)
  if (return_individual_scores == TRUE){
    return(list(combined_scores = combined_remea_scores,
                individual_db_scores = remea_indv_dbs))
  } else {
    return(combined_remea_scores)
  }
}

#' Complete ReMEA analysis of all perturbation types
#'
#' @description
#' Performs full perturbation analysis of given proteomics data and selected tumour type.
#'
#' @details
#' Returns a `list` with 8 slots. Includes ReMEA analysis for all three perturbation
#' types and averaged analysis of RNAi and CRISPR. Both combined scores and individual
#' db results are returned.
#'
#' @param protein_data `data.frame` Proteomics differential testing results. Protein name
#' must be included as protein_acc. Fold change data must be in column 2.
#' @param tumour_type `string` Whether to use markers from correlation analysis of solid tumours
#' or non-solid tumours only. If multiple passed, these will be analysed in looped format.
#' @param signature_version `string` Signature version to use.
#' @param save_xlxs `logical` Whether or not to save scores to xlxs file.
#'
#' @return ggplot object
#' @export
#'
#' @examples
complete_ReMEA_analysis <- function(protein_data,
                                    tumour_type,
                                    signature_version = "CellLines",
                                    save_xlxs = FALSE){
  # Check if object is data.table.
  if (data.table::is.data.table(protein_data)) {
    protein_data <- as.data.frame(protein_data)
  }
  # ReMEA score for drugs
  SigEnrichDrug <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                               marker_type = c("DRUG"),
                                                               tumour_type = tumour_type)
  SigEnrichDrug_combined <- ReMEA::combine_remea_scores(remea_results = SigEnrichDrug)
  message("ReMEA scoring for drugs complete.")
  # ReMEA scores for RNAi
  SigEnrichRNAi <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                                marker_type = c("RNAi"),
                                                                tumour_type = tumour_type)
  SigEnrichRNAi_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichRNAi)
  message("ReMEA scoring for RNAi complete.")
  # ReMEA scores for CRISPR
  SigEnrichCRISPR <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                                 marker_type = c("CRISPR"),
                                                                 tumour_type = tumour_type)
  SigEnrichCRISPR_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichCRISPR)
  message("ReMEA scoring for CRISPR complete.")
  # ReMEA scores for CRISPR & RNAi
  SigEnrichGENE <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                                 marker_type = c("CRISPR", "RNAi"),
                                                                 tumour_type = tumour_type)
  SigEnrichGENE_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichGENE)
  message("ReMEA scoring for RNAi/CRISPR complete.")
  scores <- list(remea_scores_drug_combined = SigEnrichDrug_combined,
                 remea_scores_drug_individual_dbs = SigEnrichDrug,
                 remea_scores_rnai_combined = SigEnrichRNAi_combined,
                 remea_scores_rnai_individual_dbs = SigEnrichRNAi,
                 remea_scores_crispr_combined = SigEnrichCRISPR_combined,
                 remea_scores_crispr_individual_dbs = SigEnrichCRISPR,
                 remea_scores_gene_combined = SigEnrichGENE_combined,
                 remea_scores_gene_individual_dbs = SigEnrichGENE
                 )
  if (xlxs_save) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      message("Saving ReMEA scores as xlxs file...")
      formatted_time <- format(Sys.time(), "%d-%m-%Y_%H-%M")
      openxlsx::write.xlsx(scores, paste0(formatted_time, "_ReMEA_scores.xlsx"))
      message("Excel file saved as: ", file_name)
    } else {
      warning("The 'openxlsx' package is not installed. Install it to save ReMEA results as Excel files.")
    }
  }
  return(scores)
}

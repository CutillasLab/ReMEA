#' Get averaged ReMEA scores
#'
#' @param protein_data A data.table containing proteomics differential-testing
#' results. A data.frame is also accepted and will be converted to a
#' data.table. Proteins must be annotated with their UniProt Entry Names.
#' The protein identifier column can be specified using protein_id_col; if
#' omitted, the first column is used. The numeric column to analyse, such as
#' fold change, can be specified using analysis_col; if omitted, the second
#' column is used.
#' @param tumour_type `string` Whether to use markers from correlation analysis of solid tumours
#' or non-solid tumours only. If multiple passed, these will be analysed in looped format.
#' @param marker_type `string` Perturbagen type for analysis. These are CRISPR, RNAi or DRUG.
#' In case of selecting both RNAi and CRISPR, these may be combined using the `combine_remea_score` function.
#' @param signature_version `string` Signature version to use.
#' @param protein_id_col Column for protein identities. To use package signatures,
#' these must refer to UniProt entry name of proteins. If NULL, will default to first column.
#' @param analysis_col Data column to analysis. If NULL, will default to second column.
#' @param return_individual_scores `logical` whether to return results for individual db scores.
#'
#' @return `data.table` of combined ReMEA scores if return_individual_scores == FALSE,
#' else returns list of both combined and individual scores.
#' @export
#'
#' @examples
get_ReMEA_scores <- function(protein_data,
                             tumour_type = c("pan"),
                             marker_type = "DRUG",
                             signature_version = "CellLines",
                             protein_id_col = NULL,
                             analysis_col = NULL,
                             return_individual_scores = FALSE){

  if (missingArg(protein_data)){
    stop("Please provide `protein_data` argument.")
  }

  remea_indv_dbs <- ReMEA::response_marker_enrichment_analysis(protein_data = protein_data,
                                                               marker_type = marker_type,
                                                               tumour_type = tumour_type,
                                                               protein_id_col = protein_id_col,
                                                               analysis_col = analysis_col,
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
#' @param protein_data `data.table` Proteomics differential testing results. Protein name
#' must be included as protein_acc. Fold change data must be in column 2.
#' @param tumour_type `string` Whether to use markers from correlation analysis of solid tumours
#' or non-solid tumours only. If multiple passed, these will be analysed in looped format.
#' @param signature_version `string` Signature version to use.
#' @param protein_id_col Column for protein identities. To use package signatures,
#' these must refer to UniProt entry name of proteins. If NULL, will default to first column.
#' @param analysis_col Data column to analysis. If NULL, will default to second column.
#' @param save_xlxs `logical` Whether or not to save scores to xlxs file.
#'
#' @return ggplot object
#' @export
#'
#' @examples
complete_ReMEA_analysis <- function(protein_data,
                                    tumour_type=c("pan"),
                                    signature_version = "CellLines",
                                    protein_id_col = NULL,
                                    analysis_col = NULL,
                                    save_xlxs = FALSE){
  if (missingArg(protein_data)){
    stop("Please provide `protein_data` argument.")
  }
  if (!data.table::is.data.table(protein_data)) {
    protein_data <- data.table::as.data.table(protein_data)
  }
  # ReMEA score for drugs
  SigEnrichDrug <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                               marker_type = c("DRUG"),
                                                               tumour_type = tumour_type,
                                                               protein_id_col = protein_id_col,
                                                               analysis_col = analysis_col,
                                                               signature_version = signature_version)
  SigEnrichDrug_combined <- ReMEA::combine_remea_scores(remea_results = SigEnrichDrug)
  message("ReMEA scoring for drugs complete.")
  # ReMEA scores for RNAi
  SigEnrichRNAi <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                               marker_type = c("RNAi"),
                                                               tumour_type = tumour_type,
                                                               protein_id_col = protein_id_col,
                                                               analysis_col = analysis_col,
                                                               signature_version = signature_version)
  SigEnrichRNAi_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichRNAi)
  message("ReMEA scoring for RNAi complete.")
  # ReMEA scores for CRISPR
  SigEnrichCRISPR <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                                 marker_type = c("CRISPR"),
                                                                 tumour_type = tumour_type,
                                                                 protein_id_col = protein_id_col,
                                                                 analysis_col = analysis_col,
                                                                 signature_version = signature_version)
  SigEnrichCRISPR_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichCRISPR)
  message("ReMEA scoring for CRISPR complete.")
  # ReMEA scores for CRISPR & RNAi
  SigEnrichGENE <-  ReMEA::response_marker_enrichment_analysis(protein_data = na.omit(protein_data),
                                                               marker_type = c("CRISPR", "RNAi"),
                                                               tumour_type = tumour_type,
                                                               protein_id_col = protein_id_col,
                                                               analysis_col = analysis_col,
                                                               signature_version = signature_version)

  SigEnrichGENE_combined <-  ReMEA::combine_remea_scores(remea_results = SigEnrichGENE)
  message("ReMEA scoring for RNAi/CRISPR complete.")
  scores <- list(drug_combined = SigEnrichDrug_combined,
                 drug_individual_dbs = SigEnrichDrug,
                 rnai_combined = SigEnrichRNAi_combined,
                 rnai_individual_dbs = SigEnrichRNAi,
                 crispr_combined = SigEnrichCRISPR_combined,
                 crispr_individual_dbs = SigEnrichCRISPR,
                 gene_combined = SigEnrichGENE_combined,
                 gene_individual_dbs = SigEnrichGENE
                 )
  if (save_xlxs) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      message("Saving ReMEA scores as xlxs file...")
      formatted_time <- format(Sys.time(), "%d-%m-%Y_%H-%M")
      openxlsx::write.xlsx(scores, paste0(formatted_time, "_ReMEA_scores.xlsx"))
      message("Excel file saved as: ", paste0(formatted_time, "_ReMEA_scores.xlsx"))
    } else {
      warning("The 'openxlsx' package is not installed. Install it to save ReMEA results as Excel files.")
    }
  }
  return(scores)
}

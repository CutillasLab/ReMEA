#' @title Response marker enrichment analysis
#'
#' @description Performs signture delta ranks and delta zscores iteratively across
#' selected signature DBs. Determined by pertubation type, tumour type and signature
#' set.
#' @param protein_data A data.table containing proteomics differential-testing
#' results. A data.frame is also accepted and will be converted to a
#' data.table. Proteins must be annotated with their UniProt Entry Names.
#' The protein identifier column can be specified using protein_id_col; if
#' omitted, the first column is used. The numeric column to analyse, such as
#' fold change, can be specified using analysis_col; if omitted, the second
#' column is used.
#' @param marker_type `string` Perturbagen type for analysis. These are CRISPR, RNAi or DRUG.
#' Multiple markers may be passed. In case of selecting both RNAi and CRISPR, these
#' may be combined using the `combine_remea_score` function
#' @param tumour_type `string` Tumour type of signature set. If using the "CellLines"
#' signatures, set to "pan".
#' @param signature_version `string` Signature version to use.
#' @param protein_id_col Column for protein identities. To use package signatures,
#' these must refer to UniProt entry name of proteins. If NULL, will default to first column.
#' @param analysis_col Data column to analysis. If NULL, will default to second column.
#' @param dataset_selection `logical` Whether to select specific signature DBs.
#' @param proteomics_datasets Optional character vector of proteomics dataset
#'   names or patterns used to filter selected signature databases. Matching is
#'   case-insensitive and is performed against signature database names. Only used
#'   when `dataset_selection = TRUE`.
#' @param perturbation_datasets Optional character vector of perturbation dataset
#'   names or patterns used to filter selected signature databases. Matching is
#'   case-insensitive and is performed against signature database names. Only used
#'   when `dataset_selection = TRUE`.
#' @return `data.table` ReMEA scores and p values for selected signature DBs.
#' @importFrom foreach %do%
#' @import doFuture
#' @export
#'
#' @examples
#' response_marker_enrichment_analysis(protein_data,
#' marker_type=c("DRUG", "RNAi", "CRISPR"),
#' tumour_type=c("pan"))
response_marker_enrichment_analysis <- function(protein_data,
                                                marker_type=c("DRUG", "RNAi", "CRISPR"),
                                                tumour_type=c("pan"),
                                                signature_version = "CellLines",
                                                protein_id_col = NULL,
                                                analysis_col = NULL,
                                                dataset_selection = FALSE,
                                                proteomics_datasets = NULL,
                                                perturbation_datasets = NULL){
  if (missingArg(protein_data)){
    stop("Please provide `protein_data` arguments.")
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
                                         signature.type = s,
                                         protein_id_col = protein_id_col,
                                         analysis_col = analysis_col))
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
#' @param protein_data A data.table containing proteomics differential-testing
#' results. A data.frame is also accepted and will be converted to a
#' data.table. Proteins must be annotated with their UniProt Entry Names.
#' The protein identifier column can be specified using protein_id_col; if
#' omitted, the first column is used. The numeric column to analyse, such as
#' fold change, can be specified using analysis_col; if omitted, the second
#' column is used.
#' @param signatures `data.table` Signature DB with response protein signatures.
#' This can be found in the ReMEA NAMESPACE. e.g. ReMEA::list_of_signature_lists
#' @param signature.type `string` Naming variable for signatures. Main purpose if for
#' response_marker_enrichment_analysis, where this indicates signature DB ids.
#' @param protein_id_col Column that indicates protein IDs. For package signatures, this
#' must be UniProt entry name. If NULL, will default to first column.
#' @param analysis_col Data column that is to be analysed. If NULL, will default to second column.
#'
#' @return `data.table` with enrichment results
#' @export
#'
#' @examples
#' signature_enrichemnt(protein_data = protein_data, signatures  = signature_dt, signature.type = "signature_id")
signature_enrichment <- function(protein_data,
                                 signatures,
                                 signature.type="signature_id",
                                 protein_id_col = NULL,
                                 analysis_col = NULL){

  if (missingArg(protein_data) | missingArg(signatures)){
    stop("Please provide `protein_data` and `signaures` arguments.")
  }
  .resolve_col <- function(x, dt, default_i, arg_name) {
    if (is.null(x)) {
      return(names(dt)[default_i])
    }

    if (is.numeric(x)) {
      if (length(x) != 1L || x < 1L || x > ncol(dt)) {
        stop(arg_name, " must be a valid single column position.", call. = FALSE)
      }
      return(names(dt)[as.integer(x)])
    }

    if (is.character(x)) {
      if (length(x) != 1L || !x %in% names(dt)) {
        stop(arg_name, " must be a valid single column name.", call. = FALSE)
      }
      return(x)
    }

    stop(arg_name, " must be NULL, a column name, or a column position.", call. = FALSE)
  }

  if (!data.table::is.data.table(protein_data)) {
    protein_data <- data.table::as.data.table(protein_data)
  } else {
    protein_data <- data.table::copy(protein_data)
  }

  protein_id_col <- .resolve_col(protein_id_col, protein_data, 1L, "protein_id_col")
  analysis_col  <- .resolve_col(analysis_col,  protein_data, 2L, "analysis_col")

  protein_data <- stats::na.omit(protein_data)
  protein_data[, protein_id := get(protein_id_col)]
  protein_data[, analysis_col := as.numeric(get(analysis_col))]
  protein_data[, analysis_col := as.numeric(scale(analysis_col))]
  protein_data[, protein_rank := rank(analysis_col)]
  mean.population <- mean(protein_data[["analysis_col"]])
  sd.population <- sd(protein_data[["analysis_col"]])

  dt_enrichment <- signatures[
    n.resistance.markers > 2 & n.sensitivity.markers > 2,
    {
      resistance.signature <- unlist(strsplit(proteins.resistance.markers, ";"))
      sensitivity.signature <- unlist(strsplit(proteins.sensitivity.markers, ";"))

      df.prot.res <- protein_data[protein_id %in% resistance.signature]
      df.prot.sen <- protein_data[protein_id %in% sensitivity.signature]

      if (nrow(df.prot.res) > 2L && nrow(df.prot.sen) > 2L) {

        resistance.markers <- paste(unique(df.prot.res[["protein_id"]]), collapse = ";")
        sensitivity.markers <- paste(unique(df.prot.sen[["protein_id"]]), collapse = ";")

        ks.p <- ks.test(
          df.prot.res[["analysis_col"]],
          df.prot.sen[["analysis_col"]]
        )

        pvalue.ks <- ks.p$p.value
        stat.ks <- as.numeric(ks.p$statistic)

        pvalue.bws <- tryCatch({
          BWStest::bws_test(
            df.prot.res[["analysis_col"]],
            df.prot.sen[["analysis_col"]],
            method = "BWS",
            alternative = "two.sided"
          )$p.value
        }, error = function(e) {
          1
        })

        zscore.resistance <- (mean(df.prot.res[["analysis_col"]]) - mean.population) / sd.population
        zscore.sensitivity <- (mean(df.prot.sen[["analysis_col"]]) - mean.population) / sd.population

        av.rank.resistance <- mean(df.prot.res[["protein_rank"]])
        av.rank.sensitivity <- mean(df.prot.sen[["protein_rank"]])

        geommean.rank.resistance <- ReMEA::gm_mean(df.prot.res[["protein_rank"]])
        geommean.rank.sensitivity <- ReMEA::gm_mean(df.prot.sen[["protein_rank"]])

        med.rank.resistance <- median(df.prot.res[["protein_rank"]])
        med.rank.sensitivity <- median(df.prot.sen[["protein_rank"]])

        n.resistance.markers <- nrow(df.prot.res)
        n.sensitivity.markers <- nrow(df.prot.sen)

        .(
          zscore.resistance = zscore.resistance,
          zscore.sensitivity = zscore.sensitivity,
          pvalue.ks = pvalue.ks,
          stat.ks = stat.ks,
          pvalue.bws = pvalue.bws,
          resistance.markers = resistance.markers,
          sensitivity.markers = sensitivity.markers,
          av.rank.resistance = av.rank.resistance,
          av.rank.sensitivity = av.rank.sensitivity,
          geommean.rank.resistance = geommean.rank.resistance,
          geommean.rank.sensitivity = geommean.rank.sensitivity,
          med.rank.resistance = med.rank.resistance,
          med.rank.sensitivity = med.rank.sensitivity,
          n.resistance.markers = n.resistance.markers,
          n.sensitivity.markers = n.sensitivity.markers,
          signature.type = signature.type
        )
      } else {
        NULL
      }
    },
    by = perturbagen
  ]

  return(dt_enrichment)
}


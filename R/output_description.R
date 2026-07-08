#' Print output column descriptions
#'
#' Prints a formatted table describing the columns returned by the package output
#' tables. The function can describe the combined score output, the individual
#' database score output, or both.
#'
#' @details
#' The combined score output contains one row per perturbagen and summarises
#' evidence across signature selections. This includes the average effect,
#' combined Kolmogorov-Smirnov and BWS p-values, adjusted p-values, and the
#' `alpha` score.
#'
#' The individual database score output contains one row per perturbagen and
#' signature source. It includes resistance and sensitivity z-scores, rank-based
#' summaries, marker counts, test statistics, p-values, adjusted p-values, and
#' derived delta scores.
#'
#' When `verbose = TRUE`, the descriptions are printed as formatted console
#' tables and the underlying description data frame is returned invisibly. When
#' `verbose = FALSE`, nothing is printed and the description data frame is
#' returned directly.
#'
#' @param combined_scores Logical scalar. Should descriptions for the combined
#'   score output be included? Defaults to `TRUE`.
#' @param individual_scores Logical scalar. Should descriptions for the
#'   individual database score output be included? Defaults to `TRUE`.
#' @param verbose Logical scalar. Should the description table be printed to the
#'   console? Defaults to `TRUE`. If `FALSE`, the description data frame is
#'   returned without printing.
#'
#' @return
#' A data frame with the following columns:
#' \describe{
#'   \item{score_type}{The output type described by the row, for example
#'   `"Combined scores"` or `"Individual database scores"`}.
#'   \item{column}{The name of the output column.}
#'   \item{type}{The expected R type of the column.}
#'   \item{description}{A plain-language description of the column.}
#' }
#'
#' If `verbose = TRUE`, the data frame is returned invisibly. If
#' `verbose = FALSE`, the data frame is returned directly.
#'
#' @examples
#' # Print descriptions for both output types
#' output_description()
#'
#' # Print only the combined score descriptions
#' output_description(
#'   combined_scores = TRUE,
#'   individual_scores = FALSE
#' )
#'
#' # Print only the individual database score descriptions
#' output_description(
#'   combined_scores = FALSE,
#'   individual_scores = TRUE
#' )
#'
#' # Return the descriptions as a data frame without printing
#' desc <- output_description(verbose = FALSE)
#'
#' # Return only the combined score descriptions as a data frame
#' combined_desc <- output_description(
#'   combined_scores = TRUE,
#'   individual_scores = FALSE,
#'   verbose = FALSE
#' )
#'
#' @export
output_description <- function(combined_scores = TRUE,
                               individual_scores = TRUE,
                               verbose = TRUE) {

  combined_desc <- data.frame(
    column = c(
      "perturbagen",
      "av.effect",
      "combined_ks_pvalue",
      "combined_bws_pvalue",
      "alpha",
      "ks_qvalue",
      "ks_padj_bonferroni",
      "bws_qvalue",
      "bws_padj_bonferroni"
    ),
    type = c(
      "character",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    description = c(
      "The perturbation in question.",
      "The composite score across signature selections.",
      paste(
        "The combined Kolmogorov-Smirnov p-value from the different",
        "signature selections. P-values are combined using Stouffer's method."
      ),
      paste(
        "The combined BWS p-value from the different signature selections.",
        "P-values are combined using Stouffer's method."
      ),
      "-log10(combined_ks_pvalue) * av.effect.",
      "Benjamini-Hochberg corrected combined_ks_pvalue.",
      "Bonferroni corrected combined_ks_pvalue.",
      "Benjamini-Hochberg corrected combined_bws_pvalue.",
      "Bonferroni corrected combined_bws_pvalue."
    ),
    stringsAsFactors = FALSE
  )

  individual_desc <- data.frame(
    column = c(
      "perturbagen",
      "zscore.resistance",
      "zscore.sensitivity",
      "pvalue.ks",
      "stat.ks",
      "pvalue.bws",
      "resistance.markers",
      "sensitivity.markers",
      "av.rank.resistance",
      "av.rank.sensitivity",
      "geommean.rank.resistance",
      "geommean.rank.sensitivity",
      "med.rank.resistance",
      "med.rank.sensitivity",
      "n.resistance.markers",
      "n.sensitivity.markers",
      "signature.type",
      "delta.zscore",
      "total_counts",
      "delta_rank_mean",
      "delta_rank_median",
      "delta_rank_geomean",
      "max_delta_score",
      "ks_qvalue",
      "ks_p_padj_bonf",
      "bws_qvalue",
      "bws_padj_bonf"
    ),
    type = c(
      "character",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "character",
      "character",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "integer",
      "integer",
      "character",
      "numeric",
      "integer",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    description = c(
      paste(
        "The perturbation in question. For compounds, this includes the",
        "site-specific code from the GDSC2 drug annotations."
      ),
      "The z-score for the resistance signature.",
      "The z-score for the sensitivity signature.",
      paste(
        "The p-value from the Kolmogorov-Smirnov test comparing the",
        "resistance and sensitivity signatures."
      ),
      paste(
        "The Kolmogorov-Smirnov test statistic comparing the sensitivity",
        "and resistance signatures."
      ),
      paste(
        "The BWS p-value comparing the resistance and sensitivity",
        "signatures."
      ),
      "Proteins used in the resistance signature.",
      "Proteins used in the sensitivity signature.",
      "Average rank of the resistance signature.",
      "Average rank of the sensitivity signature.",
      "Geometric mean of the rank of the resistance signature.",
      "Geometric mean of the rank of the sensitivity signature.",
      "Median rank of the resistance signature.",
      "Median rank of the sensitivity signature.",
      "Number of resistance markers.",
      "Number of sensitivity markers.",
      "Origin of the signature, structured as dataset_cancertype_modality.",
      "zscore.sensitivity - zscore.resistance.",
      "n.resistance.markers + n.sensitivity.markers.",
      "log2(av.rank.sensitivity) - log2(av.rank.resistance).",
      "log2(med.rank.sensitivity) - log2(med.rank.resistance).",
      paste(
        "log2(geommean.rank.sensitivity) -",
        "log2(geommean.rank.resistance)."
      ),
      paste(
        "Maximum absolute delta score across delta.zscore,",
        "delta_rank_geomean, delta_rank_mean, and delta_rank_median."
      ),
      paste(
        "Benjamini-Hochberg corrected Kolmogorov-Smirnov p-value.",
        "Adjusted within each signature.type."
      ),
      paste(
        "Bonferroni corrected Kolmogorov-Smirnov p-value.",
        "Adjusted within each signature.type."
      ),
      paste(
        "Benjamini-Hochberg corrected BWS p-value.",
        "Adjusted within each signature.type."
      ),
      paste(
        "Bonferroni corrected BWS p-value.",
        "Adjusted within each signature.type."
      )
    ),
    stringsAsFactors = FALSE
  )

  outputs <- list()

  if (isTRUE(combined_scores)) {
    outputs[["Combined scores"]] <- combined_desc
  }

  if (isTRUE(individual_scores)) {
    outputs[["Individual database scores"]] <- individual_desc
  }

  if (length(outputs) == 0) {
    out <- data.frame(
      score_type = character(),
      column = character(),
      type = character(),
      description = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(
      rbind,
      lapply(names(outputs), function(nm) {
        cbind(
          score_type = nm,
          outputs[[nm]],
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(out) <- NULL
  }

  if (!isTRUE(verbose)) {
    return(out)
  }

  if (nrow(out) == 0) {
    message("No output column descriptions selected.")
    return(invisible(out))
  }

  for (nm in names(outputs)) {
    cat("\n", nm, "\n", sep = "")
    .print_description_table(outputs[[nm]])
  }

  invisible(out)
}

#' Print a wrapped console description table
#'
#' Internal helper used by [output_description()] to print a data frame of column
#' descriptions as a wrapped console table.
#'
#' @param x A data frame with columns `column`, `type`, and `description`.
#' @param width Integer. Console width used to determine wrapping. Defaults to
#'   `getOption("width", 120)`.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of printing a
#'   formatted table to the console.
#'
#' @keywords internal
.print_description_table <- function(x, width = getOption("width", 120)) {
  stopifnot(is.data.frame(x))

  headers <- c(
    column = "Column",
    type = "Type",
    description = "Description"
  )

  x <- x[names(headers)]

  column_width <- max(nchar(c(headers[["column"]], x$column)), na.rm = TRUE)
  type_width <- max(nchar(c(headers[["type"]], x$type)), na.rm = TRUE)

  description_width <- width - column_width - type_width - 10
  description_width <- max(35, description_width)

  col_widths <- c(
    column = column_width,
    type = type_width,
    description = min(
      max(nchar(c(headers[["description"]], x$description)), na.rm = TRUE),
      description_width
    )
  )

  wrap_cell <- function(value, width) {
    value <- as.character(value)

    if (is.na(value) || !nzchar(value)) {
      return("")
    }

    wrapped <- strwrap(value, width = width)
    if (length(wrapped) == 0) {
      ""
    } else {
      wrapped
    }
  }

  pad <- function(value, width) {
    value <- as.character(value)
    paste0(value, strrep(" ", max(0, width - nchar(value))))
  }

  border <- paste0(
    "+-",
    paste(
      vapply(col_widths, function(w) strrep("-", w), character(1)),
      collapse = "-+-"
    ),
    "-+"
  )

  print_row <- function(values) {
    wrapped <- Map(wrap_cell, values, col_widths)
    row_height <- max(lengths(wrapped))

    for (i in seq_along(wrapped)) {
      length(wrapped[[i]]) <- row_height
      wrapped[[i]][is.na(wrapped[[i]])] <- ""
    }

    for (line in seq_len(row_height)) {
      cat(
        "| ",
        paste(
          vapply(
            names(col_widths),
            function(col) pad(wrapped[[col]][line], col_widths[[col]]),
            character(1)
          ),
          collapse = " | "
        ),
        " |\n",
        sep = ""
      )
    }
  }

  cat(border, "\n", sep = "")
  print_row(as.list(headers))
  cat(border, "\n", sep = "")

  for (i in seq_len(nrow(x))) {
    print_row(as.list(x[i, , drop = FALSE]))
  }

  cat(border, "\n", sep = "")
}

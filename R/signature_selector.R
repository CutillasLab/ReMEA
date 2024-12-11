#' Select signature
#'
#' @param signature_version `string` ID of signature list version
#'
#' @return ReMEA signature list.
#' @export
#'
#' @examples
select_signature <- function(signature_version) {
  if (signature_version %in% names(ReMEA::list_of_signature_lists)) {
    return(ReMEA::list_of_signature_lists[[signature_version]])
  }
  return(NULL)  # Return NULL if the combination is not found
}

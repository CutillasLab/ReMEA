#' Geometric mean
#'
#' @param x numeric vector.
#' @param na.rm `logical` whether or not to exclude na values.
#'
#' @return `float` Geometric mean of x.
#' @export
#'
#' @examples
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

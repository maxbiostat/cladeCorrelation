#' Logarithm of the difference between two logarithms
#'
#' @param x log of first value.
#' @param y log of second value.
#'
#' @return log of the difference of exp(x) and exp(y).
#' @export log_diff_exp
#'
#'
log_diff_exp <- function(x, y) {
  # if(x <= y) stop("computing the log of a negative number")
  if(y == -Inf){
    return (x)
  }else{
    return (x + log1p(-exp(y-x)) )
  }
}
#' Returns the size of a clade.
#'
#' @param x a string containing a clade.
#'
#' @return a numeric value with the size of x
#' @export get_clade_size
#'
#' @examples
#' all.c4 <- make_all_clades(4)
#' get_clade_size(all.c4)
get_clade_size <- function(x){
  length(get_clade_elements(x))
}
get_clade_size <- Vectorize(get_clade_size)

#' Get elements of a clade
#'
#' @param x a string containing a clade.
#'
#' @return a vector of strings with the elements of \code{x}.
#' @export get_clade_elements
#'
#' @examples
#' all.c4 <- make_all_clades(4)
#' get_clade_elements(all.c4)
get_clade_elements <- function(x){
  y <- gsub("\\{", "", gsub("\\}", "", x))
  strsplit(y, ",")[[1]]
}
get_clade_elements <- Vectorize(get_clade_elements)

#' Returns whether two clades are compatible
#'
#' @param x a clade.
#' @param y a clade.
#'
#' @return logical. \code{TRUE} if \code{x} and \code{y} are compatible and \code{FALSE} otherwise.
#' @export compatible
#'
#' @examples
#' compatible(x = "{t1,t2}", y = "{t1,t3,t4}")
#' compatible(x = "{t1,t2}", y = "{t1,t2,t4}")
compatible <- function(x, y){
  x.el <- get_clade_elements(x)
  y.el <- get_clade_elements(y)
  sx <- length(x.el)
  sy <- length(y.el)
  Inters <- intersect(x.el, y.el)
  invCond <- length(Inters) > 0 && length(Inters) != min(sx, sy)
  return(!invCond)
}
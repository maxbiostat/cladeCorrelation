#' Computes the joint probability for two clades for a given number of taxa.
#'
#' @param c1 a string containing a clade.
#' @param c2 a string containing a clade.
#' @param n number of taxa.
#' @param log logical. If \code{TRUE} computes the log probability (default is \code{FALSE}).
#' @param verbose logical. If \code{TRUE} it prints the case (default is \code{FALSE}).
#'
#' @return (log) joint probability.
#' @export joint_prob_clades
#'
#' @examples
#' joint_prob_clades(c1 = "{t1,t2}", c2 = "{t1,t2}", n = 5) ## case 1
#' joint_prob_clades(c1 = "{t1,t2}", c2 = "{t1,t2,t4}", n = 5) ## case 2/3
#' joint_prob_clades(c1 = "{t1,t2}", c2 = "{t3,t4,t5}", n = 5) ## case 4
#' joint_prob_clades(c1 = "{t1,t2}", c2 = "{t3,t4}", n = 5) ## case 5
#' joint_prob_clades(c1 = "{t1,t2}", c2 = "{t1,t3,t4}", n = 5) ## case 6
joint_prob_clades <- function(c1, c2, n, log = FALSE, verbose = FALSE){
  c1.el <- get_clade_elements(c1)
  c2.el <- get_clade_elements(c2)
  a <- length(c1.el)
  b <- length(c2.el)
  if(identical(sort(c1.el), sort(c2.el))){
    Case <- 1
  }else{
    if(!compatible(c1, c2)){
      Case <- 6
    }else{
      Inters <- intersect(c1.el, c2.el)
      if(length(Inters) > 0){
        Case <- 2
      }else{
        if(a + b == n){
          Case <- 4
        }else{
          Case <- 5
        }
      }
    }
  }
  if(verbose)  cat("case", Case, "\n")
  return(pclades(a = a, b = b, n = n, case = Case, log = log))
}
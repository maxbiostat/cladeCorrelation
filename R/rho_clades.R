#' Computes the correlation between two clades.
#'
#' @param c1 a string containing a clade.
#' @param c2 a string containing a clade.
#' @param n number of taxa.
#'
#' @return a data.frame containing the marginal probabilities, joint probability and correlation.
#' @export rho_clades
#'
#' @examples
#' rho_clades(c1 = "{t1,t2}", c2 = "{t1,t2}", n = 5) ## case 1
#' rho_clades(c1 = "{t1,t2}", c2 = "{t1,t2,t4}", n = 5) ## case 2/3 
#' rho_clades(c1 = "{t1,t2}", c2 = "{t3,t4,t5}", n = 5) ## case 4
#' rho_clades(c1 = "{t1,t2}", c2 = "{t3,t4}", n = 5) ## case 5
#' rho_clades(c1 = "{t1,t2}", c2 = "{t1,t3,t4}", n = 5) ## case 0
rho_clades <- function(c1, c2, n){
  c1.el <- get_clade_elements(c1)
  c2.el <- get_clade_elements(c2)
  a <- length(c1.el)
  b <- length(c2.el)
  pAB <- joint_prob_clades(c1 = c1, c2 = c2, n = n)
  pA <- pclade(k = a, n = n)
  pB <- pclade(k = b, n = n)  
  covariance <- (pAB - pA*pB)
  if(any(a==1, b==1)){
    rho <- 1
  }else{
    rho <-  sign(covariance) * exp(
      log(abs(covariance))-
        0.5*(log(pA) + log1p(-pA) + log(pB) + log1p(-pB))
    )
  }
  return(data.frame(pi = pA, pj = pB, pij = pAB, cov = covariance, rho = rho))
}

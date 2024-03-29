#' Create all (non-trivial) clades on \code{n} taxa.
#'
#' Creates all clades on \code{n} taxa in the form of a vector of  strings.
#' Clades are of the form \code{"{t1,t4}"}, \code{"{t1,t2,t4}"}.
#' @param n number of taxa.
#'
#' @return a vector of strings containing all clades on \code{n} taxa.
#' @importFrom utils combn
#' @export make_all_clades
#' @details  Excludes 'trivial' clades such as \code{"{t1}"} or \code{"{t1,t2,t3,t4,t5}"} (for 5 taxa).
#' @examples
#' make_all_clades(4)
make_all_clades <- function(n){
  tips <- paste("t", 1:n,  sep = "")
  inner <- 2:(n-1)
  Pos <- lapply(inner, function(k) combn(n, k))
  AllClades <- lapply(Pos,
                      function(cn) {
                        apply(cn, 2, function(x) paste("{", paste("t", x, collapse = ",", sep = "") , "}",sep = ""))
                      } ) 
  ans <- unlist(AllClades)
  if(length(ans) != 2^n-n-2) stop("Did not generate the correct number of clades")
  return(ans)
}

#' Construct a grid of 'all vs all' correlations for a list of K clades on \code{n} taxa.
#'
#' @param Clades a character vector with all clades to be included.
#' @param n number of taxa. 
#' @param diagonal logical. If \code{TRUE} includes the diagonal entries,
#'  which should all be 1 (default is \code{FALSE}).
#' @param ncores number of cores to be used.
#'
#' @return a grid of correlations between all clades.
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @export make_clade_corr_grid
#'
#' @examples
#' clades.n4 <- make_all_clades(4)
#' make_clade_corr_grid(clades.n4, n = 4)
make_clade_corr_grid <- function(Clades, n, diagonal = FALSE, ncores = 2){
  Clades <- gsub(" ", "", Clades) # avoid complications by removing spaces
  Var1 <- Var2 <- NULL
  K <- length(Clades)
  if(diagonal){
    posGrid <- subset(expand.grid(Var1 = 1:K, Var2 = 1:K), Var1 <= Var2) 
  }else{
    posGrid <- subset(expand.grid(Var1 = 1:K, Var2 = 1:K), Var1 < Var2)
  }
  Rho <- parallel::mclapply(1:nrow(posGrid),
                            function(i) rho_clades(c1 = Clades[[posGrid[i, 1]]],
                                                   c2 = Clades[[posGrid[i, 2]]], n = n),
                            mc.cores = ncores)
  names(posGrid) <- c("i", "j")
  return(cbind(posGrid, data.table::rbindlist(Rho)))
}

#' Build the matrix of joint probabilities for a list of K clades on \code{n} taxa.
#'
#' @param Clades a character vector with all K clades to be included.
#' @param n number of taxa.
#' @param ncores number of cores to be used.
#'
#' @return a  K x  K matrix of joint probabilities.
#' @export make_clade_joint_mat
#'
#' @examples
#' clades.n4 <- make_all_clades(4)
#' make_clade_joint_mat(Clades = clades.n4, n = 4)
make_clade_joint_mat <- function(Clades, n, ncores = 2){
  Grid <- make_clade_corr_grid(Clades = Clades,
                               n = n,
                               diagonal = TRUE,
                               ncores = ncores)
  K <- length(Clades)
  M <- matrix(0, ncol = K, nrow = K)
  for (k in 1:nrow(Grid)){
    M[Grid[k, ]$i, Grid[k, ]$j] <- Grid[k, ]$pij
    M[Grid[k, ]$j, Grid[k, ]$i] <- Grid[k, ]$pij
  } 
  return(M)
}

#' Build the variance-covariance matrix for a list of K clades on \code{n} taxa.
#' 
#' @param Clades a character vector with all K clades to be included.
#' @param n number of taxa
#' @param ncores number of cores to be used.
#'
#' @return a K x K matrix of covariances.
#' @export make_clade_cov_mat
#'
#' @examples
#' clades.n4 <- make_all_clades(4)
#' make_clade_cov_mat(Clades = clades.n4, n = 4)
make_clade_cov_mat <- function(Clades, n, ncores = 2){
  Grid <- make_clade_corr_grid(Clades = Clades,
                               n = n,
                               diagonal = TRUE,
                               ncores = ncores)
  K <- length(Clades)
  M <- matrix(0, ncol = K, nrow = K)
  for(k in 1:nrow(Grid)){
    M [Grid[k, ]$i, Grid[k, ]$j] <- Grid[k, ]$cov
    M [Grid[k, ]$j, Grid[k, ]$i] <- Grid[k, ]$cov
  } 
  return(M)
}

#' Build the correlation matrix for a list of K clades on \code{n} taxa.
#'
#' @param Clades a character vector with all clades to be included.
#' @param n number of taxa.
#' @param ncores number of cores to be used.
#'
#' @return a  K x  K correlation matrix.
#' @export  make_clade_corr_mat
#'
#' @examples
#' clades.n4 <- make_all_clades(4)
#' make_clade_corr_mat(Clades = clades.n4, n = 4)
make_clade_corr_mat <- function(Clades, n, ncores = 2){
  Grid <- make_clade_corr_grid(Clades = Clades,
                               n = n,
                               ncores = ncores)
  K <- length(Clades)
  M <- matrix(0, ncol = K, nrow = K)
  for(k in 1:nrow(Grid)) M [Grid[k, ]$i, Grid[k, ]$j] <- Grid[k, ]$rho
  M <- M + t(M)
  diag(M) <- 1
  return(M)
}

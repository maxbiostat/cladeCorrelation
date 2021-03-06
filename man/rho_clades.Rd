% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rho_clades.R
\name{rho_clades}
\alias{rho_clades}
\title{Computes the correlation between two clades.}
\usage{
rho_clades(c1, c2, n)
}
\arguments{
\item{c1}{a string containing a clade.}

\item{c2}{a string containing a clade.}

\item{n}{number of taxa.}
}
\value{
a data.frame containing the marginal probabilities, joint probability and correlation.
}
\description{
Computes the correlation between two clades.
}
\examples{
rho_clades(c1 = "{t1,t2}", c2 = "{t1,t2}", n = 5) ## case 1
rho_clades(c1 = "{t1,t2}", c2 = "{t1,t2,t4}", n = 5) ## case 2/3 
rho_clades(c1 = "{t1,t2}", c2 = "{t3,t4,t5}", n = 5) ## case 4
rho_clades(c1 = "{t1,t2}", c2 = "{t3,t4}", n = 5) ## case 5
rho_clades(c1 = "{t1,t2}", c2 = "{t1,t3,t4}", n = 5) ## case 0
}

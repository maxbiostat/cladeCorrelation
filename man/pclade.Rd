% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclade.R
\name{pclade}
\alias{pclade}
\title{Probability of a clade of size \code{k} on \code{n} taxa.}
\usage{
pclade(k, n, log = FALSE)
}
\arguments{
\item{k}{size of clade.}

\item{n}{number of taxa.}

\item{log}{logical if TRUE (default is FALSE) returns log probability.}
}
\value{
(log) probability of a clade of size \code{k} for \code{n} taxa.
}
\description{
Probability of a clade of size \code{k} on \code{n} taxa.
}
\details{
Lemma 4.2 of Zhu, Degnan & Steel (2011).
}
\examples{
pclade(k = 2, n = 5)
}
\references{
Zhu, S., Degnan, J. H., & Steel, M. (2011).
Clades, clans, and reciprocal monophyly under neutral evolutionary models.
Theoretical population biology, 79(4), 220-227.
}

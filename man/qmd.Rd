% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{qmd}
\alias{qmd}
\title{qmd() compute the quadratic mean diameter of a stand or plot.}
\usage{
qmd(dbh, expansion)
}
\arguments{
\item{dbh}{: double | vector of diameter at breast height}

\item{expansion}{: double | vector of expansion factors}
}
\value{
Returns the qmd.
}
\description{
Compute the quadratic mean diameter (basal area weighted mean diameter) of a stand or plot.
}
\examples{
data(treelist)
# compute the quadratic mean diameter
qmd( treelist$dbh, treelist$tpa )

}

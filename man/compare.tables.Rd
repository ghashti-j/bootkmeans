\name{compare.tables}
\alias{compare.tables}
\title{Contingency tables comparing true labels to fitted clusterings}
\description{
Given the output of \code{\link{compare.clusters}} and a vector of true class labels,
prints confusion tables for: (i) hard kmeans labels, (ii) the ootstrap augmented kmeans
MAP out-of-bag labels, and (optionally) (iii) fuzzy cmeans hard labels.
}

\usage{
compare.tables(full.res = NULL, true.labs = NULL)
}

\arguments{
  \item{full.res}{A list returned by \code{\link{compare.clusters}}, containing components
  \code{km}, \code{bkm}, and \code{fkm} (\code{fkm} on if argument \code{what = "all"} in function
  \code{compare.clusters}.}
  \item{true.labs}{A vector of true class labels.}
}

\details{
For kmeans, hard labels are taken from \code{full.res$km$cluster}. For bootstrap kmeans, labels
are taken from \code{full.res$bkm$clusters}. If \code{full.res$what == "all"} results are also
taken from \code{full.res$fkm$clus}, which are the hard cluster assignmented from the fuzzy cmeans
algorithm.

The function prints two or three contingency tables to the console, with three presented if
\code{\link{compare.clusters}} has argument \code{what = "all"}, and two otherwise.
}

\author{
 Jeffrey L. Andrews \email{jeff.andrews@ubc.ca} and Jesse S. Ghashti \email{jesse.ghashti@ubc.ca}
}

\seealso{
\code{\link{boot.kmeans}}, \code{\link{compare.clusters}}, \code{\link{bootk.hardsoftvis}}, \code{\link[stats]{kmeans}}, \code{\link[fclust]{FKM}}
}

\examples{
set.seed(1)
x <- as.matrix(iris[, -5])

# fit three methods (kmeans, bootstrap kmeans, fuzzy cmeans)
res <- compare.clusters(x, groups = 3, nstart = 10, what = "all")

# compare contigency tables
compare.tables(res, true.labs = iris$Species)
}

\keyword{clustering}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{bootstrap}
\keyword{kmeans}
\keyword{fuzzy}



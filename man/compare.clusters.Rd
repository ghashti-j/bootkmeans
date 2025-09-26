\name{compare.clusters}
\alias{compare.clusters}

\title{Compare Traditional kmeans, bootstrap augmented kmeans, and fuzzy cmeans}

\description{
Fits three clustering procedures on the same data: standard \code{\link[stats]{kmeans}},
our bootstrap augmented kmeans algorithm \code{\link{boot.kmeans}}, and (optionally) fuzzy cmeans
from \code{\link[fclust]{FKM}}. Returns the fitted objects of all three whose object can be passed into
\code{\link{compare.clusters}} to compare side-by-side confusion matrices.
}

\usage{
compare.clusters( data = NULL,
                  groups = NULL,
                  seed = 13462,
                  nstart = 50,
                  what = "all")
}

\arguments{
  \item{data}{Numeric matrix or data frame of row observations and column variables. Required.}
  \item{groups}{Number of clusters \eqn{K}. Required.}
  \item{seed}{Optional integer random seed for reproducibility.}
  \item{nstart}{Number of random starts for initialization for all methods.}
  \item{what}{Character flag; if \code{"all"} (default), include fuzzy cmeans (\code{fkm}) in the output.}
}

\details{
The function runs the following algorithms:
\itemize{
  \item \code{km}: \code{stats::kmeans(data, centers = groups, nstart = nstart)}.
  \item \code{bkm}: \code{boot.kmeans(data, groups, nstart = nstart, returnall = FALSE)}.
  \item \code{fkm} (if \code{what == "all"}): \code{fclust::FKM(data, k = groups, RS = nstart)}.
}
}

\value{
A named list with components:
\item{km}{\code{kmeans} fit object.}
\item{bkm}{\code{"BSKMeans"} object returned by \code{boot.kmeans}.}
\item{fkm}{(Only if \code{what == "all"}) \code{fclust} fuzzy cmeans fit.}
\item{what}{Echo of the \code{what} argument.}
}

\references{
Ghashti, J.S., Andrews, J.L., Thompson, J.R.J., Epp, J. and H.S. Kochar (2025). A bootstrap augmented $k$-means algorithm for fuzzy partitions. Submitted.

Bezdek, J.C. (1981). \emph{Pattern recognition with fuzzy objective function algorithms}. New York: Plenum.

Hartigan, J.A. and M.A. Wong (1979). Algorithm AS 136: A K-means clustering algorithm. \emph{Applied Statistics, 28}, 100–108.

Ferraro, M.B., Giordani P. and A. Serafini (2019). fclust: An R Package for Fuzzy Clustering, \emph{The R Journal, 11}.
}

\author{
 Jeffrey L. Andrews \email{jeff.andrews@ubc.ca} and Jesse S. Ghashti \email{jesse.ghashti@ubc.ca}
}

\seealso{
\code{\link{boot.kmeans}}, \code{\link{compare.tables}}, \code{\link{bootk.hardsoftvis}}, \code{\link[stats]{kmeans}}, \code{\link[fclust]{FKM}}
}

\examples{
set.seed(1)
x <- as.matrix(iris[, -5])

# compare all three methods
res <- compare.clusters(x, groups = 3, nstart = 10, what = "all")

# hard clusters from bootstrap kmeans
table(res$bkm$clusters, iris$Species)

# fuzzy memberships from fuzzy cmeans
head(res$fkm$U)

# compare class labels
cbind(res$bkm$clusters[1:5], res$fkm$clus[1:5,2], res$km$cluster[1:5])
}

\keyword{clustering}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{bootstrap}
\keyword{kmeans}
\keyword{fuzzy}

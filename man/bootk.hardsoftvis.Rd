\name{bootk.hardsoftvis}
\alias{bootk.hardsoftvis}
\title{Visualize hard vs. soft assignments from bootstrap kmeans}

\description{
Plots the results of \code{\link{boot.kmeans}} highlighting which observations are assigned with full certainty (hard) versus fractional out-of-bag membership (soft/fuzzy). Either produces a full scatterplot matrix using all variables or a 2D scatterplot of chosen variables.
}

\usage{
bootk.hardsoftvis(data = NULL, res, plotallvars = FALSE, var1 = NULL, var2 = NULL)
}

\arguments{
  \item{data}{Numeric data frame or matrix used for clustering in \code{\link{boot.kmeans}}. Required.}
  \item{res}{Result list returned from \code{\link{boot.kmeans}} (an object of class \code{"BSKMeans"}).}
  \item{plotallvars}{Logical; if \code{TRUE}, plot all pairwise scatterplots via \code{\link[graphics]{pairs}},
  otherwise \code{FALSE} requires \code{var1} and \code{var2} arguments for a 2D scatterplot. Default \code{FALSE}.}
  \item{var1}{Integer column number for the x-axis variable when \code{plotallvars = FALSE}.}
  \item{var2}{Integer column number for the y-axis variable when \code{plotallvars = FALSE}.}
}

\details{
Each observation is classified as \emph{hard} if any entry of its membership row \code{U[i,]} is exactly 1, and \emph{soft} otherwise.
These categories are mapped to colors green for hard assignments, blue for soft/fuzzy. With \code{plotallvars = TRUE}, a scatterplot matrix of all variables is drawn.
With \code{plotallvars = FALSE}, only the two specified variables are plotted, with axis labels taken from the column names of \code{data}.
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

# run bootstrap kmeans
res <- boot.kmeans(data = x, groups = 3, iterations = 20)

# scatterplot matrix of all variables
bootk.hardsoftvis(x, res, TRUE)

# scatterplot matrix of variable 1 and variable 2
bootk.hardsoftvis(x, res, plotallvars = FALSE, var1 = 1, var2 = 2)
}

\keyword{clustering}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{bootstrap}
\keyword{kmeans}
\keyword{fuzzy}


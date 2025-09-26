\name{boot.kmeans}
\alias{boot.kmeans}

\title{
Bootstrap augmented kmeans algorithm for fuzzy partitions
}

\description{
Repeatedly bootstraps the rows of a data matrix, runs \code{\link[stats]{kmeans}} on each resample (with optional seeding for given centres), tracks per-observation allocations using squared Euclidean distance, and aggregates results into out-of-bag (OOB) fuzzy memberships, hard clusters, and averaged cluster centres. Iterations can stop adaptively using a serial-correlation test on the objective trace.
}

\usage{
boot.kmeans(
  data = NULL,
  groups = NULL,
  iterations = 500,
  nstart = 1,
  export = FALSE,
  display = FALSE,
  pval = 0.05,
  itermax = 10,
  maxsamp = 1000,
  verbose = FALSE,
  returnall = FALSE
)
}

\arguments{
  \item{data}{Numeric matrix or data frame of row observations and column variables. Required.}
  \item{groups}{Either and integer number of clusters \eqn{K}; or a \eqn{K \times p} numeric matrix of initial centres. Required.}
  \item{iterations}{Initial number of bootstrap iterations to run before considering stopping (\code{default = 500}).}
  \item{nstart}{Passed to \code{\link[stats]{kmeans}} when \code{groups} is an integer (number of random starts, \code{default = 1}).}
  \item{export}{Logical; if \code{TRUE}, saves a JPEG of the objective trace at each iteration (\code{plot<i>.jpg}). Defaults to \code{FALSE}.}
  \item{display}{Logical; if \code{TRUE}, plots the most recent objective values during fitting. Defaults to \code{FALSE}.}
  \item{pval}{Significance threshold for adaptive stopping. When the Breusch–Godfrey test p-value on the last \code{iterations} objective values is not below \code{pval}, the procedure stops.}
  \item{itermax}{Maximum number of iterations per kmeans run (passed to \code{kmeans(iter.max = ...)}).}
  \item{maxsamp}{Upper bound on total iterations if adaptive stopping keeps extending (\code{default = 1000}).}
  \item{verbose}{Logical; if \code{TRUE}, print iteration counter and latest test p-value while running. Defaults to \code{FALSE}.}
  \item{returnall}{Logical; if \code{TRUE}, return full per-iteration objects (centres, kmeans fits, OOB lists); otherwise a smaller object of final results if return. Defauls to \code{TRUE}.}
}

\details{
Each iteration draws a bootstrap sample of rows, runs \code{kmeans} on the resample (first using either supplied centres or \code{nstart} random starts; subsequent iterations use the previous iteration's centres),
and computes squared Euclidean distances from every original observation to each current centre using \code{\link[stats]{mahalanobis}} with the identity
covariance. Observations are allocated to their nearest centre and these allocations are tracked across iterations.

Out-of-bag (OOB) sets are the observations note included in a given bootstrap sample. For each observation, its OOB allocations across
the most recent \code{iterations} runs are tallied to produce a fuzzy membership matrix (\eqn{U}) and a hard label by maximum membership.

Convergence is assessed adaptively: on the trace of summed per-observation minimum squared distances (the kmeans objective) over the most recent
\code{iterations} runs, a Breusch–Godfrey serial-correlation test (\code{\link[lmtest]{bgtest}} applied to a regression of the objective on
iteration index) is computed. If the p-value is below \code{pval} and \code{iterations < maxsamp}, one more iteration is added; otherwise the
loop terminates. Final centres are the elementwise mean of the centres over the last \code{iterations} runs.

}

\value{
An object of class \code{"BSKMeans"}: a list with components
\item{U}{\eqn{n \times K} matrix of OOB fuzzy cluster memberships.}
\item{clusters}{Integer vector of length \eqn{n} of hard cluster labels.}
\item{centres}{\eqn{K \times p} matrix of averaged centres over the last \code{iterations} runs.}
\item{p.value}{Final Breusch–Godfrey test p-value used for stopping.}
\item{iterations}{Total number of iterations actually run.}
\item{occurences}{\eqn{n \times \text{iterations}} matrix of per-iteration allocations for all observations.}
\item{size}{Number of clusters \eqn{K}.}
\item{soslist}{Numeric vector of objective values by iteration.}
\item{centrelist}{(If \code{returnall = TRUE}) list of per-iteration centre matrices; otherwise \code{NULL}.}
\item{ooblist}{(If \code{returnall = TRUE}) list of OOB index vectors by iteration; otherwise \code{NULL}.}
\item{kmlist}{(If \code{returnall = TRUE}) list of \code{kmeans} fit objects by iteration; otherwise \code{NULL}.}
}

\references{
Ghashti, J.S., Andrews, J.L., Thompson, J.R.J., Epp, J. and H.S. Kochar (2025). A bootstrap augmented $k$-means algorithm for fuzzy partitions. Submitted.

Breusch, T.S. (1978). Testing for Autocorrelation in Dynamic Linear Models, \emph{Australian Economic Papers}, 17, 334-355.

Godfrey, L.G. (1978). Testing Against General Autoregressive and Moving Average Error Models when the Regressors Include Lagged Dependent Variables', \emph{Econometrica}, 46, 1293-1301.
}

\author{
 Jeffrey L. Andrews \email{jeff.andrews@ubc.ca} and Jesse S. Ghashti \email{jesse.ghashti@ubc.ca}
}

\seealso{
\code{\link{compare.clusters}}, \code{\link{compare.tables}}, \code{\link{bootk.hardsoftvis}}, \code{\link[stats]{kmeans}}, \code{\link[lmtest]{bgtest}}
}

\examples{
set.seed(1)

# basic usage
x <- as.matrix(iris[, -5])
fit <- boot.kmeans(data = x, groups = 3, iterations = 50, itermax = 20, verbose = TRUE)
table(fit$clusters, iris$Species)

# basic usage with initial cluster centres supplied
centres.init <- x[sample(nrow(x), 3), ]
fit2 <- boot.kmeans(data = x, groups = centres.init, iterations = 50)

# plot objective trace
plot(fit$soslist, type = "l", xlab = "Iteration", ylab = "Objective Function Value")
}

\keyword{clustering}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{bootstrap}
\keyword{kmeans}
\keyword{fuzzy}

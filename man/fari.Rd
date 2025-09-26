\name{fari}
\alias{fari}

\title{Frobenius Adjusted Rand Indice for Comparing Two Partition Matrices}

\description{
Computes fuzzy generalizations of the  Adjusted Rand Index based on Frobenius inner products of membership matrices. These measures extend the classical indices to compare both hard and fuzzy partitions.
}

\usage{
fari(a, b)
}

\arguments{
  \item{a}{An \eqn{n \times G_1} matrix of hard or fuzzy cluster memberships, where each row sums to 1.}
  \item{b}{An \eqn{n \times G_2} matrix of hard or fuzzy cluster memberships, where each row sums to 1.}
}

\value{
A single numeric value
  \item{fari}{The Frobenius Adjusted Rand index between \code{a} and \code{b}.}
}

\references{
Andrews, J.L., Browne, R. and C.D. Hvingelby (2022). On Assessments of Agreement Between Fuzzy Partitions. \emph{Journal Classification, 39}, 326–342.
}


\examples{
set.seed(1)
a <- matrix(runif(600), nrow = 200, ncol = 3)
a <- a / rowSums(a)
b <- matrix(runif(600), nrow = 200, ncol = 3)
b <- b / rowSums(b)

fari(a, b)
}

\keyword{cluster}
\keyword{fuzzy}
\keyword{ari}
\keyword{rand index}

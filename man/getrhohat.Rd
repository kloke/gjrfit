\name{getrhohat}
\alias{getrhohat}
\title{Function to estimate AR(1) parameter}
\description{
Estimates AR(1) parameter for a linear model with cluster correlated errors,
where the within block error structure is assumed to follow an AR(1).
}
\usage{
getrhohat(x, y, block, visit, eps = 0.01)
}
\arguments{
  \item{x}{N x p design matrix}
  \item{y}{N x 1 response vector}
  \item{block}{N x 1 vector denoting block/cluster membership}
  \item{visit}{N x 1 vector denoting visit within a block/cluster}
  \item{eps}{estimate is set to 1-eps if estimated rhohat >= 1}
}
\details{
Uses method described in described in Terpstra, et. al. (2000, 2001).  c.f. Section 5.6 of Hettmansperger and McKean (2011).
}
\value{
Estimated value of rho, a scalar.
}
\references{
Hettmansperger, T. P. and McKean, J. W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd Edition}, Boca Raton, FL: CRC Press.

Terpstra, J, McKean, JW, and Naranjo, JD (2000),
Highly efficient weighted Wilcoxon estimates for autoregression, \emph{Statistics}, 35, 45-80.

Terpstra, J, McKean, JW, and Naranjo, JD (2001),
Weighted Wilcoxon estimates for autoregression, \emph{Australian \& New Zealand Journal of Statistics}, 43, 399-419.
}
\author{
John Kloke <kloke@biostat.wisc.edu>
}

\seealso{ \code{\link{gjrfit.ar1}} }
\keyword{robust}

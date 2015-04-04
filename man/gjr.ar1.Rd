\name{gjrfit.ar1}
\alias{gjrfit.ar1}
\title{Generalized Joint Rankings estimation with AR(1) within cluster correlation} 
\description{
Uses a transformation-retransformation approach to estimate parameters in a linear model with cluster correlated error.
The within cluster errors is assumed to follow an AR(1) process.
}
\usage{
gjrfit.ar1(x, y, block, visit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{N x p design matrix}
  \item{y}{N x 1 response vector}
  \item{block}{N x 1 vector denoting block/cluster membership}
  \item{visit}{N x 1 vector denoting visit within a block/cluster}
}
\details{
First, the AR(1) parameter is estimated using the method of Terpstra, et. al. (2000, 2001); c.f. Section 5.6 of Hettmansperger and McKean (2011).
Second, the response vector is transformated to uncorrelated (working independence) and intermediate fitted values are obtained.
Third, the intermediate fitted values are retransformated to obtain the fitted values (on the original scale).
Fourth, model parameters are obtained by solving the least squares problem.

Transformation of the response variables to working independence utilizes the eigenvalue decomposition.
The intermediate values require regression through the origin (Dixon and McKean 1996).
Inference is based on the JR estimation (Kloke, et. al. 2009; c.f. Section 5.2 of Hettmansperger and McKean) using a sandwich estimate (see Section 8.3 of Kloke and McKean 2014).

}
\value{
  \item{x}{N x p design matrix}
  \item{y}{N x 1 response vector}
  \item{block}{N x 1 vector denoting block/cluster membership}
  \item{visit}{N x 1 vector denoting visit within a block/cluster}
}
\references{
Dixon, SL and McKean (1996). Rank-based analysis of the heteroscedastic linear model,
\emph{Journal of the American Statistical Association}, 91, 699-712.

Hettmansperger, T. P. and McKean, J. W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd Edition}, Boca Raton, FL: CRC Press.

Kloke, J. and McKean, J.W. (2014), \emph{Nonparametric Statistical Methods Using R}, Boca Raton, FL: CRC Press.

Kloke, J.D., McKean, J.W., and Rashid, M. (2009), {Rank-Based Estimation and Associated Inferences for Linear Models with Cluster Correlated Errors}, \emph{Journal of the American Statistical Association}, 104, 485:384-390.

Terpstra, J, McKean, JW, and Naranjo, JD (2000),
Highly efficient weighted Wilcoxon estimates for autoregression, \emph{Statistics}, 35, 45-80.

Terpstra, J, McKean, JW, and Naranjo, JD (2001),
Weighted Wilcoxon estimates for autoregression, \emph{Australian \& New Zealand Journal of Statistics}, 43, 399-419.
}
\author{John Kloke <kloke@biostat.wisc.edu>}
\keyword{robust}

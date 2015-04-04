\name{ramus}
\alias{ramus}
\docType{data}
\title{Ramus (jaw bone) height}
\description{ Ramus (jaw bone) height of 20 boys measured at 8, 8.5, 9, 9.5 years.  
}
\usage{data("ramus")}
\format{
  A data frame with 80 observations on the following 3 variables.
  \describe{
    \item{\code{subject}}{subject id}
    \item{\code{age}}{age of subject in years}
    \item{\code{ramusht}}{ramus height in mm}
  }
}
\source{
Grizzle, J and Allen, D (1969), Analysis of growth and dose response curves.  \emph{Biometrics}, 25, 357-382.
}
\references{
Cook, T. and DeMets, D. (2008),
\emph{Introduction to Statistical Methods for Clinical Trials},
Boca Raton, FL: CRC Press.
}
\examples{
with(ramus,interaction.plot(age,subject,ramusht,xlab='Age',ylab='Ramus Height (mm)',legend=FALSE,lty=1))
}
\keyword{datasets}

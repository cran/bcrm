\name{plot.bcrm}
\alias{plot.bcrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the estimated dose-toxicity curve
}
\description{
The estimated dose-toxicity curve using the Bayesian continuous reassessment method is plotted for the patients thus far recruited into the trial
}
\usage{
\method{plot}{bcrm}(x, \ldots)
}
\arguments{
  \item{x}{
	An object of class "bcrm", as returned by \code{\link{bcrm}}
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
The estimated 2.5\%, 25\%, 50\%, 75\%, 97.5\% quantiles of the probability of toxicity are plotted for each dose. Additionally, a histogram of the number of toxicities and non-toxicities is plotted at each experimented dose.
}
\author{
Michael Sweeting \email{michael.sweeting@mrc-bsu.cam.ac.uk}
}
\seealso{
\code{\link{bcrm}}
}

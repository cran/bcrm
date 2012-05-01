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
\method{plot}{bcrm}(x, file = NULL, each = FALSE, \ldots)
}
\arguments{
  \item{x}{
	An object of class "bcrm", as returned by \code{\link{bcrm}}
}
  \item{each}{
	Should posterior summaries be plotted after each recruited cohort? Defaults to FALSE.
}
  \item{file}{
File name where the dose-response plots are stored, in a pdf format. The program will ammend the current sample size to the end of the file name.
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

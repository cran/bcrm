\name{plot.bcrm.sim}
\alias{plot.bcrm.sim}
\title{
Plot the operating characteristics from the simulated trials
}
\description{
Plots of the operating characteristics obtained from a CRM simulation.
}
\usage{
\method{plot}{bcrm.sim}(x, trajectories = FALSE, file = NULL, \ldots)
}
\arguments{
  \item{x}{
	An object of class "bcrm.sim", as returned by \code{\link{bcrm}} when conducting a simulation.
}
  \item{trajectories}{
	Should a summary plot of the trajectories of administered dose levels be plotted? Defaults to FALSE. 
}
  \item{file}{
File name where the operating characteristic plot is stored, in a pdf format.
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
This function plots the sample size distribution (if variable), the experimentation distribution, the recommended dose distribution and the percentage of subjects who experience the toxicity outcome (dose-limiting toxicity).
If \code{trajectories = TRUE} then summary statistics of administered dose levels for each patient are plotted instead.
}
\author{
Michael Sweeting \email{michael.sweeting@mrc-bsu.cam.ac.uk}
}
\seealso{
\code{\link{print.bcrm.sim}}, \code{\link{bcrm}}
}

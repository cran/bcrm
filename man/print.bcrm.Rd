\name{print.bcrm}
\alias{print.bcrm}
\alias{print.bcrm.sim}
\title{
Print information regarding a trial conducted using the Bayesian continuous reassessment method
}
\description{
Print method for a trial or series of trials conducted using a \code{\link{bcrm}} model. 
}
\usage{
\method{print}{bcrm}(x, \ldots)
\method{print}{bcrm.sim}(x, tox.cutpoints = NULL, trajectories = FALSE, \ldots)
}
\arguments{
  \item{x}{
An object of class "bcrm" or "bcrm.sim" as returned by \code{\link{bcrm}}
}
 \item{tox.cutpoints}{
An optional argument passed to \code{print.bcrm.sim} specifying the cutpoints of toxicity for which the operating characteristics are to be categorised. Defaults to \code{seq(from=0,to=1,by=0.2)}
}
 \item{trajectories}{
Should the individual simulation dose and outcome trajectories be returned? Defaults to FALSE. 
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
If a single trial is conducted, then the \code{\link{print}} function currently produces summary information about the design used, the data observed, current posterior estimates of toxicity, and the next recommended dose level. 
If a simulation study is conducted, then the following operating characteristics are printed:
\describe{
	\item{Experimentation proportion}{Proportion of patients recruited to each dose, and to each true region of toxicity, across the simulated trials}
	\item{Recommendation proportion}{Proportion of trials that recommend each of the dose levels as the final maximum tolerated dose (i.e. with toxicity "closest" to the target toxicity level), and the associated regions of true toxicity for the recommended MTDs}
}
If \code{trajectories = TRUE} then the dose level administered and outcome observed are returned as matrices for every patient (column) in every simulation (row). 
}


\author{
Michael Sweeting \email{michael.sweeting@mrc-bsu.cam.ac.uk}
}
\seealso{
\code{\link{bcrm}}
}

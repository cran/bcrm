\name{threep3}
\alias{threep3}
\title{
Calculate all possible trial pathways for the standard 3+3 design, together with their probability of occurring
}
\description{
 All possible pathways of a standard 3+3 design are calculated and assigned a probability of occurring. This facilitates the calculation of operating characteristics, using \code{\link{print.threep3}} and \code{\link{plot.threep3}}. 
}
\usage{
threep3(truep, start = 1, dose = NULL)
}
\arguments{
  \item{truep}{
A vector of length \code{k} (the number of doses being considered in the trial), with values equal to the true probabilities of toxicity at the dose levels.
}
 \item{start}{Starting dose level. Defaults to 1, i.e. the lowest dose level}
 \item{dose}{Optional vector of length \code{k} of actual doses for presentation purposes
}
}
\details{
The dose-escalation schema used here relates to that defined by Chang et al. (2006), which incorporates dose de-escalation also. Variations of this design exist in the literature (see Storer 1989, Reiner et al. 1999).

The first cohort of three patients are administered the starting dose (usually the lowest dose). For any cohort in the trial, escalation/de-escalation proceeds as follows:
\itemize{
	\item If none of the three patients last treated experience a DLT, then dose the next three patients at the next highest dose level;
	\item If one of the three patients last treated experiences a DLT, then dose the next three patients at the current dose level;
	\item The trial stops when at least two out of three to six patients experience DLTs at a specific dose level AND six patients have been dosed at the dose level below that given to the last cohort. The MTD is deemed to be the dose level below that given to the last cohort. If less than six patients have been dosed at the dose level below, dose a further three patients at this level.
	\item If at least two patients in the first dose level experience a DLT the trial is stopped for safety;
}
If dose-escalation extends to doses outside of that defined by \code{dose}, the MTD is determined to be the largest dose in \code{dose}.

}
\value{
\code{threep3} returns an object of class "threep3".
The function \code{\link{print}} (i.e. \code{\link{print.threep3}}) can be used to obtain operating characteristics of the design used. 

An object of class "threep3" is a list with the following components:
\item{prob}{A vector with the probabilities of each design occurring. As all possible designs are calculated, this vector sums to one}
\item{ssize}{A vector with the sample size of each design}
\item{mtd}{A vector of dose levels giving the recommended maximum tolerated dose (MTD) at the end of the trial}
\item{exp}{A vector of length \code{k} giving the average trial experimentation proportions at each dose level}
\item{dlt.no}{A vector with the number of toxicities (DLTs) that occur in each trial}
\item{truep}{The true probabilities of toxicity at each dose level, specified by the user}
}

\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. \url{http://www.jstatsoft.org/v54/i13/}

 Chang A., Ganz P., Hayes D., Kinsella T., Pass H., Schiller J., Stone R., Strecher V. \emph{Oncology: An Evidence-Based Approach}. Springer (2006).

 Storer B. Design and Analysis of Phase I Clinical Trials. \emph{Biometrics} (1989) 45: 925--937.

 Reiner E., Paoletti X., O'Quigley J. Operating characteristics of the standard phase I clinical trial design. \emph{Computational Statistics & Data Analysis} (1999) 30: 303--315.

 Neuenschwander B., Branson M., Gsponer T. Critical aspects of the Bayesian approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27: 2420--2439.

}

\author{
Graham Wheeler \email{graham.wheeler@mrc-bsu.cam.ac.uk} and

Michael Sweeting \email{michael.sweeting@mrc-bsu.cam.ac.uk}
}
\seealso{
\code{\link{threep3}}
}

\examples{
## What are the operating characteristics of a standard 3+3 design if we conside only the first 
## 12 doses of the dose-escalation cancer trial example as described in Neuenschwander et al 2008.
## Pre-defined doses
dose<-c(1,2.5,5,10,15,20,25,30,40,50,75,100)
## Pre-specified probabilities of toxicity
p.tox0<-c(0.010,0.015,0.020,0.025,0.030,0.040,0.050,0.100,0.170,0.300,0.400,0.500)

\dontrun{
design.threep3<-threep3(p.tox0,dose)
print(design.threep3)
plot(design.threep3)
}
}

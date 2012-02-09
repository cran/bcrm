\name{bcrm}
\alias{bcrm}

\title{
Bayesian continuous reassessment method (CRM) designs for Phase I dose escalation trials.
}
\description{
Implements a wide variety of Bayesian CRM designs, including 1-parameter, 2-parameter and Escalation With Overdose Control (EWOC) designs. The program can run interactively, allowing the user to enter outcomes after each cohort has been recruited, or via simulation to assess operating characteristics.
}
\usage{
bcrm(N, tox = NULL, notox = NULL, p.tox0 = NULL, sdose = NULL, dose = NULL, 
     ff, prior.alpha, cohort = 3, target.tox, constrain = TRUE,
     sdose.calculate = "mean", pointest = "plugin", tox.cutpoints = NULL,
	loss = NULL, start = NULL, simulate = FALSE, nsims = 1, truep = NULL, 
	method = "exact", burnin.itr = 2000, production.itr = 2000,
     bugs.directory = "c:/Program Files/WinBUGS14/", plot = FALSE, file = NULL)
}

\arguments{
  \item{N}{
The final sample size of the trial 
}
  \item{tox}{
A vector of length \code{k} listing the number of patients who have experienced the outcome (toxicity) at each dose level \code{1,...,k}. If missing, then it is assumed that no data have thus far been collected.
}
  \item{notox}{
A vector of length \code{k} listing the number of patients who have not experienced the outcome (toxicity) at each dose level \code{1,...,k}. If missing, then it is assumed that no data have thus far been collected.
}
  \item{p.tox0}{
A vector of length \code{k} listing the prior probabilities of experiencing the outcome at each dose level \code{1,...k}. The standardised dose levels (CRM "skeleton") are formed from these probabilities using the inverse of the functional form, with a plug-in estimate for the prior mean or median of alpha, as specified in \code{ff}, \code{prior.alpha} and \code{sdose.calculate}. Alternatively standardised dose levels can be given directly using \code{sdose}.
}
  \item{sdose}{
A vector of length \code{k} listing the standardised doses to be used in the CRM model. Only required if \code{p.tox0} is missing.
}
  \item{dose}{
Optional vector of length \code{k} of actual doses for plotting purposes
}
  \item{ff}{
A string indicating the functional form of the dose-response curve. Options are
    \describe{
     \item{ht}{ 1-parameter hyperbolic tangent}
      \item{logit1}{ 1-parameter logistic}
      \item{power}{ 1-parameter power}   
	\item{logit2}{ 2-parameter logistic}
    } 
}
  \item{prior.alpha}{
A list of length 3 containing the distributional information for the prior. The first element is a number from 1-4 specifying the type of distribution. Options are
    \enumerate{
      \item Gamma(a,b), where a=shape, b=scale: mean=a*b, variance=a*b*b
      \item Uniform(a,b), where a=min, b=max
      \item Lognormal(a,b), where a=mean on the log scale, b=standard deviation on the log scale
	\item Bivariate Lognormal(a,b), where a=mean vector on the log scale, b=Variance-covariance matrix on the log scale. This prior should be used only in conjunction with a two-parameter logistic model.
    } 
The second and third elements of the list are the parameters a and b, respectively.
}
  \item{cohort}{
The size of each cohort of patients that are sequentially recruited to the trial. Defaults to 3
}
  \item{target.tox}{
The target toxicity probability. Defaults to 1/3. 
}
  \item{constrain}{
Should a dose-skipping constraint be placed on the escalation procedure, as imposed by a modified CRM? Defaults to TRUE.
}
  \item{sdose.calculate}{
What plug-in estimate of the prior alpha should be used to calculate the standardised doses? Options are \code{"mean"} (default) or \code{"median"}. Only required if \code{sdose} is missing.
}
  \item{pointest}{
Which summary estimate of the posterior distribution should be used to choose the next dose. Options are \code{"plugin"} (default) where the posterior mean of the model parameter(s) is plugged into the function form to obtain estimates of toxicity, or \code{"mean"} where the posterior mean probabilities of toxicity are directly used. Alternatively, a number between 0 and 1 can be specified representing the quantile of the maximum tolerated dose (MTD) posterior distribution (e.g. 0.5 specifies the posteior median). This produces an Escalation With Overdose Control (EWOC) design if the quantile is less than 0.5 (see details). Currently, EWOC designs must be fit using MCMC methods.
}
  \item{tox.cutpoints}{
A vector of cutpoints for toxicity intervals if these are to be used to choose next dose. Defaults to NULL. For example Underdosing [0,0.2], Target dosing (0.2, 0.35], Excessive toxicity (0.35, 0.60], Unacceptable toxicity (0.60, 1.00] set \code{tox.cutpoints=c(0.2,0.35,0.60)}. 
}
  \item{loss}{
A vector of length \code{length(tox.cutpoints)+1} specifying the losses associated with each toxicity interval, e.g. Underdosing = 1, Target dosing =0, Excessive toxicity=1, Unacceptable toxicity=2
}
  \item{start}{
Dose level used at the beginning of the trial, or if data have been provided using \code{tox} and \code{notox} this is the current dose level for the last treated cohort. Required if \code{constrain=TRUE}.
}
  \item{simulate}{
Should a simulation be conducted to assess operating characteristics? Defaults to \code{TRUE}. If \code{FALSE}, a single CRM trial is run interactively, allowing the user to input outcomes after each cohort is recruited.
}
  \item{nsims}{
Number of simulations to perform if \code{simulate==TRUE} (defaults to 1).
}
  \item{truep}{
A vector of length k giving the true probabilities of the outcome (toxicity) at each dose level \code{1,...,k} in order to simulate data. Only required if \code{simulate=TRUE}
}
  \item{method}{
Optimisation method: options are \code{"exact"} (the default), \code{"BRugs"}, or \code{"R2WinBUGS"}.
}
  \item{burnin.itr}{
Number of burn-in iterations (default 2000).
}
  \item{production.itr}{
Number of production iterations (default 2000).
}
  \item{bugs.directory}{
Directory that contains the WinBUGS executable if \code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".
}
  \item{plot}{
Should the dose-response curve be plotted after each cohort has been entered? Defaults to FALSE.
}
  \item{file}{
File name where the dose-response plots are stored, in a pdf format. The program will ammend the current sample size to the end of the file name.
}
}
\details{
	\code{bcrm} implements a Bayesian continuous reassessment method (CRM) (O'Quigley \emph{et al.}, 1990); an adaptive design in which cohorts of patients are sequentially recruited into a Phase I trial. A binary toxicity outcome is assumed (e.g. Dose Limiting Toxicity / No Dose Limiting Toxicity). The current cohort are given a dose "closest" to the specified target toxicity level, as estimated from the posterior distributions of toxicity at each dose level from the patients thus far recruited. If \code{pointest="mean"} then the posterior mean probability of toxicity is used to choose the next dose. If \code{pointest="plugin"}, however, the posterior mean of the model parameter(s) is plugged-into the functional form of the dose-toxicity model. To implement an EWOC design (Babb \emph{et al.}, 1998), \code{pointest} should be a quantile, \emph{q}, between 0 and 0.5. The posterior distribution of the MTD (the dose in which the probability of toxicity is equal to the target toxicity) is then calculated and the next patient is given dose closest to the \emph{q}th quantile of the MTD distribution. 

Alternatively, escalation can be based on intervals of toxicity from the posterior distribution using a loss function, see Neuenschwander \emph{et al.}, 2008. To implement this approach, the user should specify the cutpoints of the toxicity intervals using \code{tox.cutpoints} and the associated losses using \code{loss}.

The possible choice of dose-toxicity model can be specified using \code{ff}, and includes the 1-parameter hyberbolic tangent, logistic or power "working models", and the 2-parameter logistic model as follows:
\describe{
	\item{Hyperbolic Tangent}{\deqn{p(Tox|d^*)=\left[(tanh(d^*)+1)/2\right]^\alpha}{p(Tox|d*)=[(tanh(d*)+1)/2]^\alpha}}
	\item{Logistic (1-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(3+\alpha d^*)}{1+\exp(3+\alpha d^*)}}{p(Tox|d*)=exp(3+\alpha d*)/(1+exp(3+\alpha d*))}}
	\item{Power}{\deqn{p(Tox|d^*)={d^*}^\alpha}{p(Tox|d*)=d*^\alpha}}
  \item{Logistic (2-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(\log(\alpha_1)+\alpha_2 d^*)}{1+\exp(\log(\alpha_1)+\alpha_2 d^*)}}{ p(Tox|d*)=exp(\alpha_1+\alpha_2 d*)/(1+exp(\alpha_1+\alpha_2 d*))}}
  }
 where \eqn{\alpha>0} is the single positive-valued parameter for the 1-parameter models, and \eqn{\log(\alpha_1)}{log(\alpha_1)} and \eqn{\alpha_2>0} are the intercept and slope parameters of the 2-parameter model.
  
  The standardised doses, \eqn{d^*}{d*}, are specified by the user using \code{sdose}, or alternatively the prior probability of toxicity at each dose level is specified using \code{p.tox0}. If the latter is used, then the standardised doses are calculated using the inverse of the functional form and a plug-in estimate of the prior mean or median, as specified in \code{sdose.calculate}, as follows
  \deqn{d^* = f^{-1}(\code{p.tox0},\alpha= a)}{d* = f^{-1}(p.tox0,\alpha= a)}
where \eqn{f^{-1}} is the the inverse of the chosen functional form, and the parameter(s) of the model are set equal to \eqn{a}, either the prior mean or median of \eqn{\alpha}. 

Data that have already been accrued can be entered using the \code{tox} and \code{notox} arguments. A constrained CRM design can be implemented using \code{constrain=TRUE}, in which case dose-skipping is prohibited (i.e. the next cohort can only be dosed up to one dose level above the current cohort). If a constrained model is used then the starting dose must be specified using \code{start}. Alternatively, if data have already been accrued, then \code{start} must give the dose level of the most recent cohort.

The prior is set using \code{prior.alpha}. For example \code{prior.alpha=list(1,1,1)} specifies a Gamma prior with shape and scale parameters both equal to one (\emph{i.e.} an Exponential(1) distribution), whilst \code{prior.alpha=list(2,0,10)} specifies a Uniform(0,10) prior. 

The trial can be run interactively using \code{simulate=FALSE}, where the user enters the outcomes for each new cohort,  or as a simulation study when \code{simulate=TRUE}.

The default calculations use exact methods (\code{method="exact"}) to calculate the mean and quantiles for the posterior distributions, unless a 2-parameter model is specified, in which case MCMC methods must be used. The two choices for MCMC calculations are \code{method="BRugs"} or \code{method="R2WinBUGS"}, the former using the OpenBUGS software, whilst the latter uses WinBUGS. To implement these methods, users require either package to be installed on their system.
}
\value{
\code{bcrm} returns an object of class "bcrm" or "bcrm.sim"; the latter occuring when a simulation has been conducted (\code{simulate=TRUE}).
The function \code{\link{print}} (i.e. \code{\link{print.bcrm}} or \code{\link{print.bcrm.sim}}) can be used to obtain summary information about the design used, the data observed, current posterior estimates of toxicity, and the next recommended dose level. 

An object of class "bcrm" is a list with the following components:
\item{dose}{Range of doses}
\item{sdose}{Standardised doses}
\item{tox}{A vector of length \code{k} listing the number of patients who have experienced the outcome (toxicity) at each dose level \code{1,...,k}}
\item{notox}{A vector of length \code{k} listing the number of patients who have not experienced the outcome (toxicity) at each dose level \code{1,...,k}}
\item{alpha}{A list with components \code{alpha.mean}, the posterior mean of the parameter(s), \code{dose.mean}, the posterior mean probability of toxicity at each dose, \code{dose.sd}, the posterior standard deviation for probability of toxicity at each dose, \code{dose.quantiles}, the posterior quantiles for probability of toxicity at each dose.}
\item{ndose}{A list with components \code{ndose}, the dose level recommended for the next patient, \code{est}, the estimated probabilitieis of toxicity using the chosen metric (e.g. plugin, mean, quantile)}
\item{constrain}{Whether a constrained CRM design was used}
\item{start}{The starting dose for the latest run of the model if \code{constrain=TRUE}}
\item{target.tox}{The target toxicity level}
\item{ff}{A number from 1-4 identifying the functional form, 1 = Hyperbolic tangent, 2 = 1-parameter logistic, 3 = Power, 4 = 2-parameter logistic}
\item{method}{The calculation method used}
\item{pointest}{The summary estimate used to choose the next dose, see \code{pointest}}
\item{prior.alpha}{Information about the prior used for \eqn{alpha}, see \code{prior.alpha}}
\item{alldoses}{A vector of length \code{N} giving the dose level administered to each patient}

An object of class "bcrm.sim" is a list of length \code{nsims}. Each component is itself a list with components similar to those obtained from a "bcrm" object. The print function, \code{\link{print.bcrm.sim}} should be used to obtain operating characteristics from the simulation.

}
\references{
 O'Quigley J., Pepe M., Fisher L. Continuous reassessment method: a practical design for phase I clinical trials in cancer. \emph{Biometrics} (1990) 46: 33--48.

 Babb J., Rogatko A., Zacks S. Cancer phase I clinical trials: efficient dose escalation with overdose control. \emph{Statistics in Medicine} (1998) 17: 1103--1120.
       
 Neuenschwander B., Branson M., Gsponer T. Critical aspects of the Bayesian approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27: 2420--2439.
}

\note{
Currently, the re-parameterisation of the two-parameter model proposed by (Babb \emph{et al.}, 1998) is not implemented. Therefore, users wishing to implement an EWOC design should check whether their choice of prior for the model parameter(s) translates to a sensible prior for the MTD distribution before they implement the design. For example
\preformatted{
	prior.alpha<-list(1,1,1)
	ff<-"ht"
	target.tox<-0.2
	samples.alpha<-getprior(prior.alpha,2000)
	mtd<-find.x(ff,target.tox,alpha=samples.alpha)
	hist(mtd)
}

One-parameter models are designed as working models only, and should not be used with an escalation strategy based on intervals of the posterior probabilities of toxicity.
}

\author{
Michael Sweeting \email{michael.sweeting@mrc-bsu.cam.ac.uk} (MRC Biostatistics Unit, Cambridge, UK), drawing on code originally developed by J. Jack Lee and Nan Chen, Department of Biostatistics, the University of Texas M. D. Anderson Cancer Center
}

\seealso{
\code{\link{print.bcrm}}, \code{\link{print.bcrm.sim}}, \code{\link{plot.bcrm}}
}
\examples{
## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
## Pre-defined doses
dose<-c(1,2.5,5,10,15,20,25,30,40,50,75,100,150,200,250)
## Pre-specified probabilities of toxicity
## [dose levels 11-15 not specified in the paper, and are for illustration only]
p.tox0<-c(0.010,0.015,0.020,0.025,0.030,0.040,0.050,0.100,0.170,0.300,0.400,0.500,0.650
  ,0.800,0.900)
## Data from the first 5 cohorts of 18 patients
tox<-c(0,0,0,0,0,0,2,0,0,0,0,0,0,0,0)
notox<-c(3,4,5,4,0,0,0,0,0,0,0,0,0,0,0)
## Target toxicity level
target.tox<-0.30

## A 1-parameter power model is used, with standardised doses calculated using 
## the plug-in prior median
## Prior for alpha is lognormal with mean 0 (on log scale) 
## and standard deviation 1.34 (on log scale)
## The recommended dose for the next cohort if posterior mean is used
Power.LN.bcrm<-bcrm(N=18,tox=tox,notox=notox,p.tox0=p.tox0,dose=dose,ff="power"
  ,prior.alpha=list(3,0,1.34^2),target.tox=target.tox,constrain=FALSE
  ,sdose.calculate="median",pointest="mean")
print(Power.LN.bcrm)
plot(Power.LN.bcrm)

## Simulate 10 replicate trials of size 36 (cohort size 3) using this design 
## with constraint (i.e. no dose-skipping) and starting at lowest dose
## True probabilities of toxicity are set to pre-specified probabilities (p.tox0) 
Power.LN.bcrm.sim<-bcrm(N=36,p.tox0=p.tox0,dose=dose,ff="power"
  ,prior.alpha=list(3,0,1.34^2),target.tox=target.tox,constrain=TRUE
  ,sdose.calculate="median",pointest="mean",start=1,simulate=TRUE,nsims=10,truep=p.tox0)
print(Power.LN.bcrm.sim)

## A 2-parameter model, using priors as specified in Neuenschwander et al 2008.
## Posterior mean used to choose the next dose
## Standardised doses using reference dose, 250mg
sdose<-log(dose/250)
## Bivariate lognormal prior for two parameters
mu<-c(2.15,0.52)
Sigma<-rbind(c(0.84^2,0.134),c(0.134,0.80^2))
TwoPLogistic.mean.bcrm<-bcrm(N=18,tox=tox,notox=notox,sdose=sdose,dose=dose,ff="logit2"
  ,prior.alpha=list(4,mu,Sigma),target.tox=target.tox,constrain=FALSE
  ,pointest="mean",method="BRugs")
print(TwoPLogistic.mean.bcrm)
plot(TwoPLogistic.mean.bcrm)

## A 2-parameter model, using an EWOC design with feasibility bound (MTD quantile) 
## of 0.25 to choose the next dose
TwoPLogistic.EWOC0.25.bcrm<-bcrm(N=18,tox=tox,notox=notox,sdose=sdose,dose=dose
  ,ff="logit2",prior.alpha=list(4,mu,Sigma),target.tox=target.tox,constrain=FALSE
  ,pointest=0.25,method="BRugs")
print(TwoPLogistic.EWOC0.25.bcrm)
plot(TwoPLogistic.EWOC0.25.bcrm)

## A 2-parameter model, using a loss function based on intervals of toxicity to choose
## the next dose
## Toxicity cut-points
tox.cutpoints<-c(0.2,0.35,0.6)
## Losses associated with toxicity intervals 
## [0,0.2]=1, (0.2,0.35]=0, (0.35,0.6]=1, (0.6,1]=2
loss<-c(1,0,1,2)
TwoPLogistic.tox.intervals.bcrm<-bcrm(N=18,tox=tox,notox=notox,sdose=sdose,dose=dose
  ,ff="logit2",prior.alpha=list(4,mu,Sigma),target.tox=target.tox,constrain=FALSE
  ,tox.cutpoints=tox.cutpoints,loss=loss,method="BRugs")
print(TwoPLogistic.tox.intervals.bcrm)
plot(TwoPLogistic.tox.intervals.bcrm)
## Greater loss associated with overdosing and unacceptable toxicity
## [0,0.2]=1, (0.2,0.35]=0, (0.35,0.6]=2, (0.6,1]=4
loss2<-c(1,0,2,4)
TwoPLogistic.tox.intervals.2.bcrm<-bcrm(N=18,tox=tox,notox=notox,sdose=sdose,dose=dose
  ,ff="logit2",prior.alpha=list(4,mu,Sigma),target.tox=target.tox,constrain=FALSE
  ,tox.cutpoints=tox.cutpoints,loss=loss2,method="BRugs")
print(TwoPLogistic.tox.intervals.2.bcrm)
plot(TwoPLogistic.tox.intervals.2.bcrm)


}
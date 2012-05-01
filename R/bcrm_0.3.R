## Bayesian CRM - extending original code of J. Jack Lee and Nan Chen, Department of Biostatistics, the University of Texas M. D. Anderson Cancer Center
## Now using exact inference, R2WinBUGS and BRugs 

### MJS 24/02/12

# ----------------------------------------------------------------------
# 	bcrm. Conduct a Bayesian CRM trial simulating outcomes from a true model
# 	stop	     --> A list of stopping rules containing information on when the trial should be terminated
#			nmax = Maximum sample size
#			nmtd = Maximum number to be treated at final MTD estimate			
#			precision = Lower and Upper risks of toxicity that must contain the 95% posterior interval (2.5%ile and 97.5%ile) of the recommended dose in order to stop the trial
#			nmin = Minimum sample size required in the trial
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#	p.tox0     --> prior probabilities of response at dose levels
#     sdose      --> standardised dose levels (given if p.tox0 is missing)
#    dose		--> optional vector of dose labels (for printing and plotting purposes)
#	ff	     --> functional form of dose-response curve
#                    "ht" - hyperbolic tangent
#                    "logit1" - 1-parameter logistic
#                    "power" - Power
#				  "logit2" - Two-parameter logistic		 
#     prior.alpha--> list containing the information to generate the
#                   prior distribution for parameter (vector) alpha
#                   [[1]] - the prior distribution type:
#                       1 - gamma: mean=a*b; var=a*b*b
#                       2 - uniform: unif(a,b)
#					3 - lognormal: lnorm(a,b), where a=mean on the log scale, b=var on the log scale
#				     4 - log Multivariate normal(a,b), where a=mean vector on log scale, b=Variance-covariance matrix (log scale) (for multiparameter functional forms)
#                   [[2]] - a: first parameter (scalar/vector) of the prior distribution
#                   [[3]] - b: second parameter (scalar/vector) of the prior distribution
#     cohort     --> cohort size - default = 3
#     target.tox --> target toxicity 
# 	constrain --> TRUE (default) or FALSE
#    sdose.calculate -> What plug-in estimate of the prior alpha should be used to calculate the standardised doses? Options are "mean" (default) or "median"
# 	pointest   --> Which summary estimate of the posterior distribution should be used to choose next dose, "plugin" (default), "mean" or a numerical quantile between 0 and 1 (e.g. 0.5). If NULL then escalation based on posterior intervals (Loss function approach) is assumed
#	tox.cutpoints --> Specify the cutpoints for toxicity intervals if these are to be used to choose next dose, e.g. Underdosing [0,0.2], Target dosing (0.2, 0.35], Excessive toxicity (0.35, 0.60], Unacceptable toxicity (0.60, 1.00]
#	loss --> The losses associated with each toxicity interval, e.g. Underdosing = 1, Target dosing =0, Excessive toxicity=1, Unacceptable toxicity=2
#     start      --> Starting dose level or if data is provided this is the current dose level for the last treated cohort. Required if constrain=TRUE
#    simulate -->  Perform a simulation to assess operating characteristics (Default=TRUE). If FALSE, a single CRM trial is run interactively, allowing the user to input outcomes after each cohort is recruited?
#    nsims	--> No. of simulations to perform if simulate==T (defaults to 1)
# 	truep	     --> True probabilities of response at dose levels to simulate data. Only should be specified if simulate=TRUE
#	method      --> Optimisation method: options are "exact" (the default), "BRugs", or "R2WinBUGS"
#     burnin.itr --> No. of burn-in iterations
#     production.itr --> No. of production iterations
#    bugs.directory --> directory that contains the WinBUGS executable if R2WinBUGS is being used, defaults to "C:/Program Files/WinBUGS14/"
#	plot	     --> Should dose response curve be plotted after each cohort is entered? Defaults to FALSE
# 	file	     --> name of the file where the dose-response plots are stored, in a pdf format
#	N --> Final sample size (Deprecated). To be replaced with stop in future versions
# ----------------------------------------------------------------------
bcrm<-function(stop=list(nmax=NULL,nmtd=NULL,precision=NULL,nmin=NULL),tox=NULL,notox=NULL,p.tox0=NULL,sdose=NULL,dose=NULL,ff,prior.alpha,cohort=3,target.tox,constrain=TRUE,sdose.calculate="mean",pointest="plugin",tox.cutpoints=NULL,loss=NULL,
	start=NULL,simulate=FALSE,nsims=1,truep=NULL,
	method="exact",burnin.itr=2000,production.itr=2000,bugs.directory="c:/Program Files/WinBUGS14/",plot=FALSE,file=NULL,N){

	# Checks of argument inputs	
	if(missing(N) & is.null(stop$nmax) & is.null(stop$nmtd) & is.null(stop$precision)) stop("At least one stopping rule must be provided using the stop argument")
	if(!missing(N)){
		 stop$nmax<-N
		 warning("N is deprecated and users should now use the stop argument to specify the maximum sample size")
	}
	if(!(length(stop$precision) %in% c(0,2))) stop("stop$precision must be a vector of length two")
	if(!is.null(stop$nmax) & !is.null(stop$nmin)) {if(stop$nmin>stop$nmax) stop("stop$nmin must be less than stop$nmax")}
	if(missing(p.tox0) & missing(sdose)) stop("Either p.tox0 or sdose must be specified")
	if(!missing(p.tox0) & !missing(sdose)) stop("Only one of p.tox0 and sdose must be specified")
	if(sdose.calculate!="mean" & sdose.calculate!="median") stop("sdose.calculate must be either `mean' or `median'")
	if((is.character(pointest) & pointest!="mean" & pointest!="plugin") | is.numeric(pointest) & (pointest<0 | pointest>1)) stop("pointest must be either `plugin', `mean' or an EWOC feasibility quantile between 0 and 1")
	if(is.numeric(pointest) & method=="exact") stop("EWOC design must be fitted using MCMC methods")
	if(!is.null(tox.cutpoints) & method=="exact") stop("Escalation based on toxicity intervals must be fit using MCMC. Please specify either method='BRugs' or method='R2WinBUGS'")

	if(simulate & is.null(truep)) stop("truep must be specified if simulating data")
	if(simulate){ 
		plot<-FALSE
		results<-list()
	}
	if(!(method %in% c("exact","BRugs","R2WinBUGS"))) stop("method must be either `exact', `BRugs' or `R2WinBUGS'")
	## Check to see if ff is one of "ht","logit1","power","logit2"
	if((!ff %in% c("ht","logit1","power","logit2"))) stop("ff must be one of `ht', `logit1', `power' or `logit2'")

	if(ff=="logit2" & method=="exact") warning("Exact method slow for 2-parameter model, suggest using BRugs (MCMC)")
	if(constrain & is.null(start)) stop("A starting/current dose level must be specified using `start' if constrain==TRUE")
	if((!is.null(tox.cutpoints) & is.null(loss)) | (is.null(tox.cutpoints) & !is.null(loss))) stop("Both tox.cutpoints and loss must be specified to conduct escalation based on toxicity intervals")
	if(!is.null(tox.cutpoints) & length(loss)!=length(tox.cutpoints)+1) stop("The number of losses must be one more than the number of cutpoints")
	if(!is.null(tox.cutpoints)) pointest<-NULL
	if(!is.null(tox.cutpoints) & ff!="logit2") warning("One-parameter models are designed as working models only, and should not be used with an escalation strategy based on intervals of the posterior probabilities of toxicity")

     if(ff=="logit2" & (length(prior.alpha[[2]])<2 | length(prior.alpha[[3]])<2)) stop("second and third components of `prior.alpha' must be vectors of size 2")

	if(missing(sdose)){
		alpha.prior.plug<-if(prior.alpha[[1]]==1){
			ifelse(sdose.calculate=="mean",prior.alpha[[2]]*prior.alpha[[3]],median(getprior(prior.alpha, 10000)))
		} else if(prior.alpha[[1]]==2){
			0.5*(prior.alpha[[2]]+prior.alpha[[3]])
		} else if(prior.alpha[[1]]==3){
			ifelse(sdose.calculate=="mean",exp(prior.alpha[[2]]+prior.alpha[[3]]/2),exp(prior.alpha[[2]]))
		} else if(prior.alpha[[1]]==4){
			if(sdose.calculate=="mean"){exp(prior.alpha[[2]]+diag(prior.alpha[[3]])/2)} else {exp(prior.alpha[[2]])}
		} 		
		sdose<-find.x(ff,p.tox0,alpha=alpha.prior.plug)
	}
	# Checking that length of truep in simulation study is same length as sdose
	if(simulate & length(truep)!=length(sdose)) stop("Length of truep must be the same as length of sdose or p.tox0")	
	# Checking that length of dose (if specified) is same length as sdose
	if(length(dose)>0 & length(dose)!=length(sdose)) stop("Length of dose must be the same as length of sdose or p.tox0")

	## Allowing access to fast exact computation if method="exact" & simulate=TRUE & stopping rules do not depend on posterior quantiles
	if(method=="exact" & simulate & is.null(stop$precision)){ method<-"exact.sim" }
	## If method=="exact" and a two-parameter model is fitted, only relevant escalation posterior quantities are calculated (apart from trial end)
	if(method=="exact" & ff=="logit2"){ 
		method<-"exact.sim"
		if(plot){
			plot<-FALSE
			warning("Plot function not available for exact computation of 2-parameter model")
		}
	}
	k<-length(sdose)
	sim<-1	

	while(sim<=nsims){
		## Prior
		if(is.null(tox) & is.null(notox)){
			new.tox<- rep(0,k)
			new.notox<-rep(0,k)
			current<-start-1
			alpha<-if(sim==1) switch(method
			,BRugs=getprior(prior.alpha, 10000)
			,R2WinBUGS=getprior(prior.alpha,10000)
			,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
			,exact.sim=Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest)
			)
		} else { 
			new.tox<-tox
			new.notox<-notox
			current<-start 

			alpha<-if(sim==1) switch(method
			,BRugs=Posterior.BRugs(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
			,R2WinBUGS=Posterior.R2WinBUGS(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr,bugs.directory)
			,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
			,exact.sim=Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest)
			)
		}
		ncurrent<-sum(new.tox+new.notox)
		alldoses<-NULL
		alltox<-NULL
		found<-FALSE
		match<-1
		ndose<-if(sim==1){
			switch(method
				,BRugs=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss)
				,R2WinBUGS=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss)
				,exact=nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
				,exact.sim=nextdose.exact.sim(alpha,sdose,ff,target.tox,constrain,pointest,current)
				)
		} else {
			results[[1]]$ndose[[1]]
		}
		if(!simulate){
			results<-list(dose=dose,sdose=sdose,tox=new.tox,notox=new.notox,ndose=list(ndose),constrain=constrain,start=start,target.tox=target.tox,ff=ff,method=method,pointest=pointest,tox.cutpoints=tox.cutpoints,loss=loss,prior.alpha=prior.alpha)
			class(results)<-"bcrm"
		} else {
			results[[sim]]<-list(dose=dose,sdose=sdose,tox=new.tox,notox=new.notox,ndose=list(ndose),constrain=constrain,start=start,target.tox=target.tox,ff=ff,method=method,pointest=pointest,tox.cutpoints=tox.cutpoints,loss=loss,prior.alpha=prior.alpha,truep=truep)
			class(results)<-"bcrm.sim"
		}
		if(plot){
			plot(results,file)
		}
		stopped<-stop.check(stop,ncurrent,ndose,new.tox,new.notox,simulate)
		while(!stopped){
			current<-ndose[[1]]
			alldoses<-c(alldoses,rep(current,cohort))
			ncurrent<-ncurrent+cohort
			if(!simulate){
				interact<-crm.interactive(new.tox,new.notox,ncurrent,cohort,ndose,dose)
				if(interact$bk==TRUE){
					results$alldoses<-alldoses
					results$alltox<-alltox
					return(results)
				}
				y<-interact$y
				new.tox<-interact$tox
				new.notox<-interact$notox
			} else {	
				y<-rbinom(1,cohort,truep[ndose[[1]]]) 
				new.notox[ndose[[1]]]<-new.notox[ndose[[1]]]+(cohort-y)
				new.tox[ndose[[1]]]<-new.tox[ndose[[1]]]+y
			}
			alltox<-c(alltox,rep(c(0,1),c(cohort-y,y)))

			if(simulate & match<length(results) & !(nsims>1000 & method=="exact.sim")){
				repeat{
					m<-all(results[[match]]$alldoses[1:length(alldoses)]==alldoses & results[[match]]$alltox[1:length(alltox)]==alltox)
					if(m){
						ndose<-results[[match]]$ndose[[length(results[[sim]]$ndose)+1]]
						found<-TRUE
						break
					} else if(match!=(length(results)-1)){
						match<-match+1
					} else {
						match<-match+1	
						found<-FALSE			
						break
					}
				}
			} 
			if(!found){
				alpha<-switch(method
					,BRugs=Posterior.BRugs(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
					,R2WinBUGS=Posterior.R2WinBUGS(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr,bugs.directory)
					,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
					,exact.sim=Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest)
					)
				ndose<-switch(method
					,BRugs=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss)
					,R2WinBUGS=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss)
					,exact=nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
					,exact.sim=nextdose.exact.sim(alpha,sdose,ff,target.tox,constrain,pointest,current)
					)
			}
			if(!simulate){
				results$tox<-new.tox
				results$notox<-new.notox
				results$ndose[[length(results$ndose)+1]]<-ndose
				results$alldoses<-alldoses
				results$alltox<-alltox
			} else{
				results[[sim]]$tox<-new.tox
				results[[sim]]$notox<-new.notox
				results[[sim]]$ndose[[length(results[[sim]]$ndose)+1]]<-ndose
				results[[sim]]$alldoses<-alldoses
				results[[sim]]$alltox<-alltox
			}
			if(plot){
				plot(results,file)
			}
			stopped<-stop.check(stop,ncurrent,ndose,new.tox,new.notox,simulate)
		}
		if(simulate){
			cat(sim,"\n")
		}
	sim<-sim+1
	}
	return(results)
}

# ----------------------------------------------------------------------
# 	CRM INTERACTIVE. Conduct a CRM trial interactively allowing user to specify outcomes after each cohort has been recruited
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#	ncurrent --> Current no. of patients in trial
#	cohort   --> Cohort size
#    ndose    --> Proposed dose level for next cohort
#	dose     --> Dose labels for each dose
# ----------------------------------------------------------------------
crm.interactive<-function(tox,notox,ncurrent,cohort,ndose,dose){
     k <- length(tox)
		repeat {
			  if(is.null(dose)){
	     		  cat("\n\n RECOMMENDED DOSE LEVEL FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent, "IS:", ndose[[1]])
		            ans <- get.dose.level(k,ncurrent,cohort)     
			  } else {
			  	  cat("\n\n RECOMMENDED DOSE FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent, "IS:", dose[ndose[[1]]])
		            ans <- get.dose(dose,ncurrent,cohort)     
			  }
     	       if (ans==-2)
          	      ans <- ifelse(is.null(dose),ndose[[1]],dose[ndose[[1]]])
	            if (ans==0) {
				  cat("\n\n EXIT AND RETURN THE RESULTS SO FAR?")            
     	           cat("\n DO YOU REALLY WANT TO EXIT ? (Y/N)  ")
          	      yn <- readline()
              	 	 if (yn=="y" || yn=="Y") break
              		  else                    next
	     		 }  
			# y.j is the number of toxicities from the treatment of
		     # cohort patients at dose level ndose
     		cat("\n ENTER NUMBER OF TOXICITIES FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent,": ")
		     y <- get.answer(cohort)               
     		# give the user a last chance to modify the treatment and outcome
		     cat("\n\n\t\t ENTERED VALUES:")
			if(is.null(dose)){
	     		cat("\n DOSE LEVEL ...", ans)
			} else {
				cat("\n DOSE ...", ans)
			}
	     	cat("\n NO. OF TOXICIES ....",y)
	     	cat("\n PRESS `RETURN' IF OK, OR ANY OTHER KEY TO ENTER NEW VALUES ")
		     key <- readline()
     		if (nchar(key)==0) break
		}
          if (ans==0)  return(list(bk=TRUE))
		if(!is.null(dose)){
			ans<-which(dose==ans)
		}
		notox[ans]<-notox[ans]+(cohort-y)
		tox[ans]<-tox[ans]+y
		return(list(tox=tox,notox=notox,y=y,bk=FALSE))
}


# ----------------------------------------------------------------------
#     User inputs the number of toxicities from the treatment at a dose level:
#	 Must be less than or equal to the size of the cohort
# ----------------------------------------------------------------------
get.answer <- function(cohort) {   
    repeat {
        ans <- as.numeric(readline())
        if (is.na(ans)) next
        if (ans != floor(ans)) next
        if (ans>=0 & ans<=cohort) return(ans)
    }
}

# ----------------------------------------------------------------------
#     ask the user to input an integer number <= n.
#     the number is checked to belong to [-1,n] and also to be
#     an integer.  `ENTER' returns to the caller with no action,
#
#     n - biggest number to be accepted
# 	 ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose.level <- function( n ,ncurrent,cohort) {
    repeat {
        cat("\n\n ENTER DOSE LEVEL BETWEEN 1 AND ", n," FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent)
        cat("\n (`RETURN' TO ACCEPT RECOMMENDATION, 0 TO EXIT AND RETURN CURRENT RESULTS)  ")
        
        ans <- readline()
        if ( nchar(ans)==0 ) return( -2 )
        ans <- as.integer(ans)
        if (is.na(ans)) next
        if ( -1<=ans && ans<=n ) return( ans )
    }
}

# ----------------------------------------------------------------------
#     ask the user to input a dose from those given
#     `ENTER' returns to the caller with no action,
#
#     dose	- dose labels
# 	 ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose <- function( dose ,ncurrent,cohort) {
    repeat {
        cat("\n\n ENTER DOSE FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent)
	   cat("\n POSSIBLE CHOICES ARE ",dose)
        cat("\n (`RETURN' TO ACCEPT RECOMMENDATION, 0 TO EXIT AND RETURN CURRENT RESULTS)  ")
        
        ans <- readline()
        if ( nchar(ans)==0 ) return( -2 )
        if (is.na(ans)) next
        if ( ans %in% dose | ans==0) return( ans )
    }
}

# ----------------------------------------------------------------------
# generates a vector of values and prior distribution
#
#     prior.alpha --> list containing the information for prior
#               [[1]] - the prior distribution type:
#                     1 - gamma: mean=a*b; var=a*b*b
#                     2 - uniform: a+b*unif(0,1)
#				   3 - lognormal: lnorm(a,b), where mean=a, var=b
#				   4 - log Multivariate normal(a,b), where a=mean vector, b=Variance-covariance matrix
#               [[2]] - a: first parameter of the prior distribution
#               [[3]] - b: second parameter of the prior distribution
#    
# ----------------------------------------------------------------------
getprior <- function(prior.alpha, n) {
    type <- prior.alpha[[1]]
    a    <- prior.alpha[[2]]
    b    <- prior.alpha[[3]]    
    if ( type==1 ) {
        prior<-rgamma(n,a)*b        
    }
    else if ( type==2 ) {
        prior<-sapply(1:length(a),function(i){runif(n,a[i],b[i])})
    }
    else if (type==3) { 
    	   prior<-rlnorm(n,a,sqrt(b))
    }
    else if (type==4) {
	   log.prior<-rmvnorm(n,a,b)
	   prior<-exp(log.prior)
    }
    return (prior)
}

# ----------------------------------------------------------------------
#     returns the vector containing sampled data from winbug
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#				  "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#     burnin.itr  --> number of burn-in iterations
#     production.itr --> number of production iterations
# ----------------------------------------------------------------------
Posterior.BRugs <- function(tox, notox,sdose,ff, prior.alpha, burnin.itr, production.itr)
{    
  if(require(BRugs) & require(R2WinBUGS)){
    all.patient <- tox + notox
    datan <-all.patient[all.patient!=0]
    datas <-tox[all.patient!=0]
    datad <-sdose[all.patient!=0]
    k <- length(datan)
    if (k == 1)
    {
        datan <- c(datan, 0)
        datas <- c(datas, 0)
        datad <- c(datad, 0)
    }
    mydata <- list(N1 = k, s = datas,n = datan,d = datad, p1 = prior.alpha[[2]], p2 = prior.alpha[[3]])
    model.file<-if (prior.alpha[[1]] == 1)
    	{
        if (ff == "ht")
            HTGamma
        else if (ff == "logit1")
            LogisticGamma
        else if (ff == "power")
            PowerGamma
	   else stop("Functional form not currently available with specified prior distribution")
    }
    else if (prior.alpha[[1]] == 2)
    {
        if (ff == "ht")
            HTUnif
        else if (ff == "logit1")
            LogisticUnif
        else if (ff == "power")
            PowerUnif
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 3)
    {
        if (ff == "ht")
            HTLogNormal
        else if (ff == "logit1")
            LogisticLogNormal
        else if (ff == "power")
            PowerLogNormal
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 4)
    {
        if (ff == "logit2")
		TwoPLogisticLogNormal            
	   else stop("Functional form not currently available with specified prior distribution")
    }    
	path.model<-file.path(tempdir(), "model.file.txt")
	path.data<-file.path(tempdir(), "data.file.txt")
	write.model(model.file,path.model)
	bugsData(mydata,path.data)
	BRugsFit(path.model,path.data, inits=NULL, numChains = 2, parametersToSave="alpha",
	    	nBurnin = burnin.itr, nIter = production.itr/2, BRugsVerbose = FALSE,DIC=F)
	if(ff=="logit2"){
		t<-cbind(samplesSample("alpha[1]"),samplesSample("alpha[2]"))
	} else {   t<- samplesSample("alpha") }
    return(t)
  }
}

# ----------------------------------------------------------------------
#     returns the vector containing sampled data from winbug
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#     burnin.itr  --> number of burn-in iterations
#     production.itr --> number of production iterations
#	 bugs.directory --> directory that contains the WinBUGS executable, defaults to C:/Program Files/WinBUGS14/
# ----------------------------------------------------------------------
Posterior.R2WinBUGS <- function(tox, notox,sdose,ff, prior.alpha, burnin.itr, production.itr,bugs.directory)
{    
  if(require(R2WinBUGS)){
    all.patient <- tox + notox
    datan <-all.patient[all.patient!=0]
    datas <-tox[all.patient!=0]
    datad <-sdose[all.patient!=0]
    k <- length(datan)
    if (k == 1)
    {
        datan <- c(datan, 0)
        datas <- c(datas, 0)
        datad <- c(datad, 0)
    }
    mydata <- list(N1 = k, s = datas,n = datan,d = datad, p1 = prior.alpha[[2]], p2 = prior.alpha[[3]])
    ## initdat<-list(list(alpha = 0.5),list(alpha=1))
    parameters<-"alpha"
    model.file<-if (prior.alpha[[1]] == 1)
    	{
        if (ff == "ht")
            HTGamma
        else if (ff == "logit1")
            LogisticGamma
        else if (ff == "power")
            PowerGamma
	   else stop("Functional form not currently available with specified prior distribution")
    }
    else if (prior.alpha[[1]] == 2)
    {
        if (ff == "ht")
            HTUnif
        else if (ff == "logit1")
            LogisticUnif
        else if (ff == "power")
            PowerUnif
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 3)
    {
        if (ff == "ht")
            HTLogNormal
        else if (ff == "logit1")
            LogisticLogNormal
        else if (ff == "power")
            PowerLogNormal
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 4)
    {
        if (ff == "logit2")
		TwoPLogisticLogNormal            
	   else stop("Functional form not currently available with specified prior distribution")
    }    
    res<-bugs(mydata, inits=NULL, parameters, model.file,n.chains=2,n.thin=1,n.iter=burnin.itr+production.itr/2,n.burnin=burnin.itr,DIC=F,bugs.directory=bugs.directory)  
    if(is.null(dim(res$summary))){
	    res$summary<-t(as.matrix(res$summary))
	}
    if(any(res$summary[,"Rhat"]>1.01)){ 
	warning("Convergence may not have been achieved: Trace plots shown")
	par(mfrow=c(dim(res$sims.array)[[3]],1))
	for(i in 1:dim(res$sims.array)[[3]]){
		plot(res$sims.array[,1,i],type="l")
		lines(res$sims.array[,2,i],col=2)
    	}
    }
	t<-apply(res$sims.array,3,rbind)
    return(t)
  }
}

# ----------------------------------------------------------------------
#     returns the posterior mean of alpha and actualised values of toxic probabilities at each dose level
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
# ----------------------------------------------------------------------
Posterior.exact<-function(tox,notox,sdose,ff,prior.alpha){    
    all.patient <- tox + notox
    data.tox <-tox[all.patient!=0]
    data.notox <-notox[all.patient!=0]
    data.dose <-sdose[all.patient!=0]
   
    ## Following code fixes bug if no data available (prior only)
    if(length(data.dose)==0){
		data.dose<-sdose[1]
		data.tox<-data.notox<-0
    }	
    wmodel<-which.f(ff)

    prior<-switch(prior.alpha[[1]]
			,"1"=function(alpha,prior.alpha){dgamma(alpha,shape=prior.alpha[[2]],scale=prior.alpha[[3]])}
		 	,"2"=function(alpha,prior.alpha){dunif(alpha,min=prior.alpha[[2]],max=prior.alpha[[3]])}
			,"3"=function(alpha,prior.alpha){dlnorm(alpha,meanlog=prior.alpha[[2]],sdlog=sqrt(prior.alpha[[3]]))}
			,"4"=function(alpha,prior.alpha){1/(alpha[,1]*alpha[,2])*dmvnorm(log(alpha),mean=prior.alpha[[2]],sigma=prior.alpha[[3]])})

	if(prior.alpha[[1]]!=4){
		## Scaling factor to prevent likelihood getting too small
		C<-1/prod(sapply(1:length(data.dose),function(d){wmodel(data.dose[d],1)^data.tox[d]*(1-wmodel(data.dose[d],1))^data.notox[d]}))
		lik<-function(dose,tox,notox,alpha,C){
			l<-rep(1,length(alpha))
			for(d in 1:length(dose)){
				l<-l*wmodel(dose[d],alpha)^tox[d]*(1-wmodel(dose[d],alpha))^notox[d]
			}
			C*l
		}
		int.norm.constant<-function(alpha,dose,tox,notox,prior.alpha){lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		int.alpha.mean<-function(alpha,dose,tox,notox,prior.alpha){alpha*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		int.dose.mean<-function(alpha,new.dose,dose,tox,notox,prior.alpha){wmodel(new.dose,alpha)*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		int.dose.sd<-function(alpha,new.dose,dose.mean,dose,tox,notox,prior.alpha){(wmodel(new.dose,alpha)-dose.mean)^2*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}

		norm.constant<-integrate(int.norm.constant,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]
		alpha.mean<-integrate(int.alpha.mean,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
		dose.mean<-sapply(sdose,function(d){integrate(int.dose.mean,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),new.dose=d,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
		dose.sd<-sapply(1:length(sdose),function(d){sqrt(integrate(int.dose.sd,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),new.dose=sdose[d],dose.mean=dose.mean[d],dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant)})

		cdf<-function(par){integrate(int.norm.constant,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),par,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant}
		fn<-function(par,q){abs(cdf(par)-q)}
		alpha.quantiles<-sapply(c(0.025,0.25,0.5,0.75,0.975),function(q){	max.x<-seq(0,10*alpha.mean,length=100)[which.min(sapply(seq(0,10*alpha.mean,length=100),function(i){fn(i,q)}))+1]
		optimize(fn,c(0,max.x),q=q, tol = .Machine$double.eps^0.50)$minimum
		})
		dose.quantiles<-sapply(sdose,function(d){wmodel(d,sort(alpha.quantiles,TRUE))})
		rownames(dose.quantiles)<-c("2.5%","25%","50%","75%","97.5%")
	} else {
		## VECTORISED FORM OF LIK
		lik<-function(dose,tox,notox,alpha){
			l<-rep(1,nrow(alpha))
			for(d in 1:length(dose)){
				l<-l*wmodel(dose[d],alpha)^tox[d]*(1-wmodel(dose[d],alpha))^notox[d]
			}
			l
		}
		joint.posterior<-function(alpha1,alpha2,dose,tox,notox,prior.alpha){lik(dose,tox,notox,cbind(alpha1,alpha2))*prior(cbind(alpha1,alpha2),prior.alpha)}
		joint.posterior.2<-function(alpha2,alpha1,dose,tox,notox,prior.alpha){lik(dose,tox,notox,cbind(alpha1,alpha2))*prior(cbind(alpha1,alpha2),prior.alpha)}
		marginal.alpha1<-function(alpha1,dose,tox,notox,prior.alpha){sapply(alpha1,function(a1){integrate(joint.posterior.2,0,Inf,alpha1=a1,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}	
		marginal.alpha2<-function(alpha2,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(joint.posterior,0,Inf,alpha2=a2,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}
		int.alpha1.mean<-function(alpha1,dose,tox,notox,prior.alpha){alpha1*marginal.alpha1(alpha1,dose,tox,notox,prior.alpha)}
		int.alpha2.mean<-function(alpha2,dose,tox,notox,prior.alpha){alpha2*marginal.alpha2(alpha2,dose,tox,notox,prior.alpha)}
		int.dose.mean<-function(alpha1,alpha2,new.dose,dose,tox,notox,prior.alpha){wmodel(new.dose,cbind(alpha1,alpha2))*joint.posterior(alpha1,alpha2,dose,tox,notox,prior.alpha)}
		int.dose.mean.2<-function(alpha2,new.dose,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(int.dose.mean,0,Inf,alpha2=a2,new.dose=new.dose,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}
		int.dose.sd<-function(alpha1,alpha2,new.dose,dose.mean,dose,tox,notox,prior.alpha){(wmodel(new.dose,cbind(alpha1,alpha2))-dose.mean)^2*joint.posterior(alpha1,alpha2,dose,tox,notox,prior.alpha)}
		int.dose.sd.2<-function(alpha2,new.dose,dose.mean,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(int.dose.sd,0,Inf,alpha2=a2,new.dose=new.dose,dose.mean=dose.mean,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}

		norm.constant<-integrate(marginal.alpha2,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]
		alpha1.mean<-integrate(int.alpha1.mean,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
		alpha2.mean<-integrate(int.alpha2.mean,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
		alpha.mean<-c(alpha1.mean,alpha2.mean)
		dose.mean<-sapply(sdose,function(d){integrate(int.dose.mean.2,0,Inf,new.dose=d,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
		dose.sd<-sapply(1:length(sdose),function(d){sqrt(integrate(int.dose.sd.2,0,Inf,new.dose=sdose[d],dose.mean=dose.mean[d],dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant)})

		dose.quantiles<-NULL

		# Using cuhre
		#prior<-function(alpha,prior.alpha){1/(alpha[1]*alpha[2])*dmvnorm(log(alpha),mean=prior.alpha[[2]],sigma=prior.alpha[[3]])}
		#wmodel<-function(dose,alpha) {1/(1+exp(-log(alpha[1])-alpha[2]*dose))}
		## Scaling factor to prevent likelihood getting too small
		#C<-1/prod(wmodel(data.dose,c(1,1))^data.tox*(1-wmodel(data.dose,c(1,1)))^data.notox)
		#lik<-function(dose,tox,notox,alpha,C){
		#	l<-prod(wmodel(dose,alpha)^tox*(1-wmodel(dose,alpha))^notox)
		#	l
		#}
		#int.norm.constant<-function(alpha,dose,tox,notox,prior.alpha){lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		#norm.constant<-cuhre(2,1,int.norm.constant,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,lower=c(0,0),upper=c(Inf,Inf),flags=list(verbose=0))$value
		#int.dose.mean<-function(alpha,new.dose,dose,tox,notox,prior.alpha){wmodel(new.dose,alpha)*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		#dose.mean<-cuhre(2,length(sdose),int.dose.mean,new.dose=sdose,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,lower=c(0,0),upper=c(500,500))$value/norm.constant

	}
	return(list(alpha.mean=alpha.mean,dose.mean=dose.mean,dose.sd=dose.sd,dose.quantiles=dose.quantiles))
}


# ----------------------------------------------------------------------
#     Cut-down version of Posterior.exact for use with simulations and two-parameter models
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
# 	pointest   --> Which summary estimate of the posterior distribution should be used to choose next dose, "plugin" (default), "mean" or a numerical quantile between 0 and 1 (e.g. 0.5). If NULL then escalation based on posterior intervals (Loss function approach) is assumed
# ----------------------------------------------------------------------
Posterior.exact.sim<-function(tox,notox,sdose,ff,prior.alpha,pointest){    

    all.patient <- tox + notox
    data.tox <-tox[all.patient!=0]
    data.notox <-notox[all.patient!=0]
    data.dose <-sdose[all.patient!=0]
   
    ## Following code fixes bug if no data available (prior only)
    if(length(data.dose)==0){
		data.dose<-sdose[1]
		data.tox<-data.notox<-0
    }	
    wmodel<-which.f(ff)

    prior<-switch(prior.alpha[[1]]
			,"1"=function(alpha,prior.alpha){dgamma(alpha,shape=prior.alpha[[2]],scale=prior.alpha[[3]])}
		 	,"2"=function(alpha,prior.alpha){dunif(alpha,min=prior.alpha[[2]],max=prior.alpha[[3]])}
			,"3"=function(alpha,prior.alpha){dlnorm(alpha,meanlog=prior.alpha[[2]],sdlog=sqrt(prior.alpha[[3]]))}
			,"4"=function(alpha,prior.alpha){1/(alpha[,1]*alpha[,2])*dmvnorm(log(alpha),mean=prior.alpha[[2]],sigma=prior.alpha[[3]])})

	if(prior.alpha[[1]]!=4){
		## Scaling factor to prevent likelihood getting too small
		C<-1/prod(sapply(1:length(data.dose),function(d){wmodel(data.dose[d],1)^data.tox[d]*(1-wmodel(data.dose[d],1))^data.notox[d]}))
		lik<-function(dose,tox,notox,alpha,C){
			l<-rep(1,length(alpha))
			for(d in 1:length(dose)){
				l<-l*wmodel(dose[d],alpha)^tox[d]*(1-wmodel(dose[d],alpha))^notox[d]
			}
			C*l
		}
		int.norm.constant<-function(alpha,dose,tox,notox,prior.alpha){lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
		norm.constant<-integrate(int.norm.constant,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]
		if(pointest=="plugin"){
			int.alpha.mean<-function(alpha,dose,tox,notox,prior.alpha){alpha*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
			alpha.mean<-integrate(int.alpha.mean,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
			dose.mean<-NULL
		} else if(pointest=="mean"){
			int.dose.mean<-function(alpha,new.dose,dose,tox,notox,prior.alpha){wmodel(new.dose,alpha)*lik(dose,tox,notox,alpha,C)*prior(alpha,prior.alpha)}
			alpha.mean<-NULL
			dose.mean<-sapply(sdose,function(d){integrate(int.dose.mean,ifelse(prior.alpha[[1]]==2,prior.alpha[[2]],0),ifelse(prior.alpha[[1]]==2,prior.alpha[[3]],Inf),new.dose=d,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
		} 
	} else {
		## VECTORISED FORM OF LIK
		lik<-function(dose,tox,notox,alpha){
			l<-rep(1,nrow(alpha))
			for(d in 1:length(dose)){
				l<-l*wmodel(dose[d],alpha)^tox[d]*(1-wmodel(dose[d],alpha))^notox[d]
			}
			l
		}
		joint.posterior<-function(alpha1,alpha2,dose,tox,notox,prior.alpha){lik(dose,tox,notox,cbind(alpha1,alpha2))*prior(cbind(alpha1,alpha2),prior.alpha)}
		joint.posterior.2<-function(alpha2,alpha1,dose,tox,notox,prior.alpha){lik(dose,tox,notox,cbind(alpha1,alpha2))*prior(cbind(alpha1,alpha2),prior.alpha)}
		marginal.alpha1<-function(alpha1,dose,tox,notox,prior.alpha){sapply(alpha1,function(a1){integrate(joint.posterior.2,0,Inf,alpha1=a1,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}	
		marginal.alpha2<-function(alpha2,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(joint.posterior,0,Inf,alpha2=a2,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}
		int.alpha1.mean<-function(alpha1,dose,tox,notox,prior.alpha){alpha1*marginal.alpha1(alpha1,dose,tox,notox,prior.alpha)}
		int.alpha2.mean<-function(alpha2,dose,tox,notox,prior.alpha){alpha2*marginal.alpha2(alpha2,dose,tox,notox,prior.alpha)}
		int.dose.mean<-function(alpha1,alpha2,new.dose,dose,tox,notox,prior.alpha){wmodel(new.dose,cbind(alpha1,alpha2))*joint.posterior(alpha1,alpha2,dose,tox,notox,prior.alpha)}
		int.dose.mean.2<-function(alpha2,new.dose,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(int.dose.mean,0,Inf,alpha2=a2,new.dose=new.dose,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}
		int.dose.sd<-function(alpha1,alpha2,new.dose,dose.mean,dose,tox,notox,prior.alpha){(wmodel(new.dose,cbind(alpha1,alpha2))-dose.mean)^2*joint.posterior(alpha1,alpha2,dose,tox,notox,prior.alpha)}
		int.dose.sd.2<-function(alpha2,new.dose,dose.mean,dose,tox,notox,prior.alpha){sapply(alpha2,function(a2){integrate(int.dose.sd,0,Inf,alpha2=a2,new.dose=new.dose,dose.mean=dose.mean,dose=dose,tox=tox,notox=notox,prior.alpha=prior.alpha)$value})}

		norm.constant<-integrate(marginal.alpha2,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]
		if(pointest=="plugin"){
			alpha1.mean<-integrate(int.alpha1.mean,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
			alpha2.mean<-integrate(int.alpha2.mean,0,Inf,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
			alpha.mean<-c(alpha1.mean,alpha2.mean)
			dose.mean<-NULL
		} else if(pointest=="mean"){
			alpha.mean<-NULL
			dose.mean<-sapply(sdose,function(d){integrate(int.dose.mean.2,0,Inf,new.dose=d,dose=data.dose,tox=data.tox,notox=data.notox,prior.alpha=prior.alpha,rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
		}
	}
	return(list(alpha.mean=alpha.mean,dose.mean=dose.mean))
}


# ----------------------------------------------------------------------
#     find the new dose and posteria mean of alpha
#
#     samples.alpha  --> sampled data of variable(s) alpha
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#		         "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#	constrain	  --> TRUE or FALSE
#	pointest	  --> "plugin", "mean" or a quantile
#	current	  --> current dose level (of last treated cohort) - 0 if this is first cohort
#
# ----------------------------------------------------------------------
nextdose<-function(samples.alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss){
	f<-which.f(ff)
	k<-length(sdose)
	samples.sdose<-sapply(sdose,function(x){f(x,samples.alpha)})
	mean<-apply(samples.sdose,2,mean)
	median<-apply(samples.sdose,2,median)
	sd<-apply(samples.sdose,2,sd)
	quantiles<-apply(samples.sdose,2,quantile,c(0.025,0.25,0.50,0.75,0.975))

	if(!is.null(loss)){
		probs<-sapply(1:k,function(i){hist(samples.sdose[,i], c(0,tox.cutpoints,1), plot = FALSE)$counts/dim(samples.sdose)[1]})
		est<-apply(loss*probs,2,sum)
		target<-0
	} else {

		if(pointest=="plugin"){
			mean.alpha<-apply(as.matrix(samples.alpha),2,mean)
			if(length(mean.alpha)>1) mean.alpha<-t(as.matrix(mean.alpha))
		} 
		if(is.numeric(pointest)){
			mtd<-find.x(ff,ptox=target.tox,alpha=samples.alpha)
			target<-quantile(mtd,pointest)
		} else {
			target<-target.tox
		}
		est<-if(pointest=="plugin"){ f(sdose,mean.alpha)
		} else {if(pointest=="mean"){apply(samples.sdose,2,mean)
			} else { sdose }}
		}
	# Next dose	
	ndose<-if(!constrain){which.min(abs(est-target))
		} else {which.min(abs(est[1:min(current+1,k)]-target))}
	
	return(list(ndose=ndose,est=est,mean=mean,sd=sd,quantiles=quantiles,target=target))
}


# ----------------------------------------------------------------------
#     find the new dose and posteria mean of alpha
#
#     alpha     --> posterior mean of alpha and posterior mean p(DLT) at each dose level
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#	constrain	  --> TRUE or FALSE
#	pointest	  --> "mean" or "plugin"
#	current	  --> current dose level (of last treated cohort) - 0 if this is first cohort
#
# ----------------------------------------------------------------------
nextdose.exact<-function(alpha,sdose,ff,target.tox,constrain,pointest,current){
	f<-which.f(ff)
	k<-length(sdose)
	est<-if(pointest=="mean") alpha$dose.mean
		 else if(pointest=="plugin") f(sdose,alpha$alpha.mean)
		 else stop("Quantile estimation not available for exact computation, please use pointest='mean' or 'plugin'")
	mean<-alpha$dose.mean
	sd<-alpha$dose.sd
	quantiles<-alpha$dose.quantiles
	ndose<-if(!constrain){which.min(abs(est-target.tox))
		} else {which.min(abs(est[1:min(current+1,k)]-target.tox))}
	return(list(ndose=ndose,est=est,mean=mean,sd=sd,quantiles=quantiles))
}

# ----------------------------------------------------------------------
#     Cut-down version of nextdose.exact for use with simulations
#
#     alpha     --> posterior mean of alpha and posterior mean p(DLT) at each dose level
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#	constrain	  --> TRUE or FALSE
#	pointest	  --> "mean" or "plugin"
#	current	  --> current dose level (of last treated cohort) - 0 if this is first cohort
#
# ----------------------------------------------------------------------
nextdose.exact.sim<-function(alpha,sdose,ff,target.tox,constrain,pointest,current){
	f<-which.f(ff)
	k<-length(sdose)
	est<-if(pointest=="mean") alpha$dose.mean
		 else if(pointest=="plugin") f(sdose,alpha$alpha.mean)
		 else stop("Quantile estimation not available for exact computation, please use pointest='mean' or 'plugin'")
	ndose<-if(!constrain){which.min(abs(est-target.tox))
		} else {which.min(abs(est[1:min(current+1,k)]-target.tox))}
	return(list(ndose=ndose,est=est))
}


# -----------------------------------------------------------------------
#     Return the calculation models
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#			   "logit2" - Two-parameter logistic
# -----------------------------------------------------------------------
which.f <- function( ff ) {
    return( switch(ff,
                    ht=function(dose,alpha) {((tanh(dose)+1)/2)^alpha},
                    logit1=function(dose,alpha) {1/(1+exp(-3-alpha*dose))},
                    power=function(dose,alpha) {dose^alpha},
                    logit2=function(dose,alpha) {1/(1+exp(-log(alpha[,1])-alpha[,2]*dose))} ))
}

# ----------------------------------------------------------------------
#     Find the dose corresponding to a certain p(DLT)
#
#     ff         --> functional form for the probability density
#     ptox       --> desired probability of DLT
#     alpha      --> the parameter(s) of the model
# ----------------------------------------------------------------------
find.x <- function( ff, ptox, alpha ) {
    if ( ff == "ht" ) {
	  x<-atanh(2*ptox^(1/alpha)-1)
    }
    else if ( ff == "logit1" )
        x <- (qlogis(ptox)-3)/alpha
    else if ( ff == "power")
        x <- exp(log(ptox)/alpha)
	else if (ff=="logit2"){
	   if(is.vector(alpha))	alpha<-matrix(alpha,ncol=2)
	   x <- (qlogis(ptox)-log(alpha[,1]))/alpha[,2]
	}
	return( x )
}

# ----------------------------------------------------------------------
#     Check whether we have reached the stopping criteria
#
#     stop         --> list of stopping criteria (nmax, nmtd, precision, nmin)
#     ncurrent     --> current sample size
#     ndose	   --> next dose information, including posterior quantities
#     new.tox	   --> No. toxicities to date
#     new.notox	   --> No. non-toxicities to date
#     simulate     --> Is this a simulation study?
# ----------------------------------------------------------------------
stop.check<-function(stop,ncurrent,ndose,new.tox,new.notox,simulate){
	answer1<-answer2<-answer3<-FALSE
	answer4<-TRUE
	if(!is.null(stop$nmin)){
		answer4<-(ncurrent>=stop$nmin)
	} 
	if(!is.null(stop$nmax)){
		answer1<-(ncurrent>=stop$nmax)
		if(answer1 & answer4 & !simulate) cat("\n Stopping: Reached maximum sample size\n")
	}	
	if(!is.null(stop$nmtd)){
		answer2<-(new.tox+new.notox)[ndose$ndose]>=stop$nmtd
		if(answer2 & answer4 & !simulate) cat("\n Stopping: Reached maximum number at MTD estimate\n")
	}
	if(!is.null(stop$precision)){
		answer3<- stop$precision[1] <= ndose$quantiles["2.5%",ndose$ndose] & ndose$quantiles["97.5%",ndose$ndose]<= stop$precision[2]
		if(answer3 & answer4 & !simulate) cat("\n Stopping: Reached required precision for MTD estimate\n")
	}

	(answer1 | answer2 | answer3) & answer4 
}

#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm
# -----------------------------------
plot.bcrm<-function(x,file=NULL,each=FALSE,...){
	dose<-if(is.null(x$dose)) x$sdose else x$dose
	dose.label<-if(is.null(x$dose)) "Standardised dose" else "Dose"
	f<-which.f(x$ff)
	if(!each){
		if(x$method %in% c("exact","exact.sim") & x$ff=="logit2"){
			df<-data.frame(dose=dose,target.tox=x$target.tox,est=x$ndose[[length(x$ndose)]]$est)
		} else {
			df<-data.frame(dose=dose,target.tox=x$target.tox,est=x$ndose[[length(x$ndose)]]$est,mean=x$ndose[[length(x$ndose)]]$mean,q2.5=x$ndose[[length(x$ndose)]]$quantiles["2.5%",],q25=x$ndose[[length(x$ndose)]]$quantiles["25%",],q50=x$ndose[[length(x$ndose)]]$quantiles["50%",],q75=x$ndose[[length(x$ndose)]]$quantiles["75%",],q97.5=x$ndose[[length(x$ndose)]]$quantiles["97.5%",])
		}
		df2<-data.frame(dose=factor(c(rep(dose,x$tox),rep(dose,x$notox)),levels=dose),Outcome=factor(c(rep("DLT",sum(x$tox)),rep("No DLT",sum(x$notox))),levels=c("DLT","No DLT")))
	} else {
		if(x$method %in% c("exact","exact.sim") & x$ff=="logit2"){
			df<-data.frame(dose=rep(dose,length(x$ndose)),target.tox=x$target.tox,cohort=rep(0:(length(x$ndose)-1),each=length(dose)),est=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$est})),
				ndose=rep(sapply(1:length(x$ndose),function(i){dose[x$ndose[[i]]$ndose]}),each=length(dose)))
		} else {
			df<-data.frame(dose=rep(dose,length(x$ndose)),target.tox=x$target.tox,cohort=rep(0:(length(x$ndose)-1),each=length(dose)),est=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$est})),mean=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$mean})),q2.5=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$quantiles["2.5%",]})),q25=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$quantiles["25%",]})),
				q50=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$quantiles["50%",]})),q75=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$quantiles["75%",]})),q97.5=c(sapply(1:length(x$ndose),function(i){x$ndose[[i]]$quantiles["97.5%",]})),
				ndose=rep(sapply(1:length(x$ndose),function(i){dose[x$ndose[[i]]$ndose]}),each=length(dose)))
		}
		df2<-data.frame()
	}
	a<-if(is.null(x$loss)){
		if(!each){ 	
			if(x$method %in% c("exact","exact.sim") & x$ff=="logit2"){
				ggplot()+geom_point(aes(x=dose,y=est),data=df)+
				geom_hline(aes(yintercept=target.tox),data=df,col=4,linetype=2)	+xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+opts(title="Posterior point estimates \n Diamond shows next recommended dose")+
				geom_point(aes(x=dose,y=est),data=df[x$ndose[[length(x$ndose)]][[1]],],size=4,col=4,shape=9)
			} else {
				ggplot()+geom_errorbar(aes(x=dose,ymin=q2.5,ymax=q97.5),colour="red",data=df)+geom_pointrange(aes(x=dose,y=q50,ymin=q25,ymax=q75),data=df,fill="red")+
				geom_hline(aes(yintercept=target.tox),data=df,col=4,linetype=2)	+xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+opts(title="Posterior p(DLT) quantiles: 2.5%, 25%, 50%, 75%, 97.5% \n Diamond shows next recommended dose")+
				geom_point(aes(x=dose,y=q50),data=df[x$ndose[[length(x$ndose)]][[1]],],size=4,col=4,shape=9)
			}
		} else {
			if(x$method %in% c("exact","exact.sim") & x$ff=="logit2"){
				ggplot()+geom_point(aes(x=dose,y=est),data=df)+
				geom_hline(aes(yintercept=target.tox),data=df,col=4,linetype=2)	+xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+opts(title="Posterior point estimates \n Diamond shows next recommended dose")+
				geom_point(aes(x=ndose,y=q50),data=df[df$dose==df$ndose,],size=4,col=4,shape=9)+
				facet_wrap(~ cohort) 
			} else {
				ggplot()+geom_errorbar(aes(x=dose,ymin=q2.5,ymax=q97.5),colour="red",data=df)+geom_pointrange(aes(x=dose,y=q50,ymin=q25,ymax=q75),data=df,fill="red")+
				geom_hline(aes(yintercept=target.tox),data=df,col=4,linetype=2)	+xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+opts(title="Posterior p(DLT) quantiles: 2.5%, 25%, 50%, 75%, 97.5% \n Diamond shows next recommended dose")+
				geom_point(aes(x=ndose,y=q50),data=df[df$dose==df$ndose,],size=4,col=4,shape=9)+
				facet_wrap(~ cohort) 
			}
		}
	} else { 
		df.intervals<-data.frame(cohort=rep(0:(length(x$ndose)-1),each=length(x$loss)),xmin=min(dose),xmax=max(dose),ymin=c(0,tox.cutpoints),ymax=c(x$tox.cutpoints,1),Loss=x$loss)
		if(!each){
			ggplot()+geom_errorbar(aes(x=dose,ymin=q2.5,ymax=q97.5),colour="red",data=df)+geom_pointrange(aes(x=dose,y=q50,ymin=q25,ymax=q75),data=df,fill="red")+
			geom_point(aes(x=dose,y=q50),data=df[x$ndose[[length(x$ndose)]][[1]],],size=4,col=4,shape=9)+
			geom_rect(aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=Loss),data=df.intervals,alpha=0.3)+scale_fill_gradient(breaks=sort(unique(df.intervals$Loss)),high="red",low="#99ccff")+
			xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+
			opts(title="Posterior p(DLT) quantiles: 2.5%, 25%, 50%, 75%, 97.5% \n Diamond shows next recommended dose")
		} else {
			ggplot()+geom_errorbar(aes(x=dose,ymin=q2.5,ymax=q97.5),colour="red",data=df)+geom_pointrange(aes(x=dose,y=q50,ymin=q25,ymax=q75),data=df,fill="red")+
			geom_point(aes(x=ndose,y=q50),data=df[df$dose==df$ndose,],size=4,col=4,shape=9)+
			geom_rect(aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=Loss),data=df.intervals,alpha=0.3)+xlab(dose.label)+ylab("Probability of DLT")+ylim(0,1)+
			opts(title="Posterior p(DLT) quantiles: 2.5%, 25%, 50%, 75%, 97.5% \n Diamond shows next recommended dose")+
			facet_wrap(~ cohort)
		}
	}
	b<-if(nrow(df2)!=0) {
		ggplot()+geom_bar(aes(x=dose,fill=Outcome),data=df2)+xlab(dose.label)+ylab("Number")+scale_fill_hue(limits=c("DLT","No DLT"))
	} else { NULL	}

	if(!is.null(file))
		pdf(paste(file,sum(x$tox+x$notox),".pdf",sep=""),...)
	if(!each){
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(2,1)))
		vplayout<-function(x,y)	viewport(layout.pos.row=x,layout.pos.col=y)
		print(a,vp=vplayout(1,1))	
		if(!is.null(b)) print(b,vp=vplayout(2,1))	
	} else {
		print(a)
	} 
	if(!is.null(file))
		dev.off()
}

#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm.sim
# -----------------------------------
plot.bcrm.sim<-function(x,trajectories=FALSE,file=NULL,...){
	dose<-if(is.null(x[[1]]$dose)) x[[1]]$sdose else x[[1]]$dose
	dose.label<-if(is.null(x[[1]]$dose)) "Standardised dose" else "Dose"

	if(trajectories){
		## sample size
		n<-sapply(x,function(i){length(i$alldoses)})
		traj.mat<-matrix(NA,nrow=length(x),ncol=max(n))
		for(i in 1:length(x)){
			traj.mat[i,1:n[i]]<-x[[i]]$alldoses
		}
		traj.df<-data.frame(patient=rep(1:max(n),each=5),Statistic=factor(rep(c("Minimum","Lower Quartile","Median","Upper Quartile","Maximum"),max(n)),levels=c("Minimum","Lower Quartile","Median","Upper Quartile","Maximum")),traj=c(apply(traj.mat,2,quantile,na.rm=T)))
		cols<- c("Median" = "black","Lower Quartile" = "blue","Upper Quartile" = "blue", "Minimum" = "red","Maximum"="red")
		lt<-c("Median" = 1,"Lower Quartile" = 2,"Upper Quartile" = 2, "Minimum" = 4,"Maximum"=4)
		ggplot()+geom_step(aes(x=patient,y=traj,group=Statistic,colour=Statistic,linetype=Statistic),data=traj.df)+
			scale_colour_manual(values=cols)+scale_linetype_manual(values=lt)+
			xlab("Patient")+ylab("Dose Level")
		if(!is.null(file))
			ggsave(paste(file,".pdf",sep=""),...)
	} else {
		# sample size
		n<-sapply(x,function(i){length(i$alldoses)})
		df.n<-data.frame(n)
		a<-if(min(n)!=max(n)){
			ggplot()+geom_histogram(aes(x=n,y=100*..density..),data=df.n,binwidth=1)+xlab("Sample size")+ylab("Percent")+opts(title="Sample size")
		} else {NULL}

		# experimentation
		exp<-rep(dose,apply(sapply(x,function(i){(i$tox+i$notox)}),1,sum))
		df.exp<-data.frame(exp)
		b<-ggplot()+geom_histogram(aes(x=exp,y=100*..density..),data=df.exp,binwidth=1)+xlab(dose.label)+ylab("Percent")+opts(title="Experimentation")

		# recommendation
		rec<-dose[sapply(x,function(i){i$ndose[[length(i$ndose)]]$ndose})]
		df.rec<-data.frame(rec)
		c<-ggplot()+geom_histogram(aes(x=rec,y=100*..density..),data=df.rec,binwidth=1)+xlab(dose.label)+ylab("Percent")+opts(title="Recommendation")

		# observed DLTs
		obs<-sapply(x,function(i){100*sum(i$tox)/sum(i$tox+i$notox)})
		df.obs<-data.frame(obs)
		d<-ggplot()+geom_histogram(aes(x=obs,y=100*..density..),data=df.obs,binwidth=1)+xlab("Percentage of subjects with DLTs")+ylab("Percent")+opts(title="DLTs")

		if(!is.null(file))
			pdf(paste(file,".pdf",sep=""),...)
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(2,2)))
		vplayout<-function(x,y)	viewport(layout.pos.row=x,layout.pos.col=y)
		if(!is.null(a)) print(a,vp=vplayout(1,1))	
		print(b,vp=vplayout(1,2))
		print(c,vp=vplayout(2,1))	
		print(d,vp=vplayout(2,2))	
		if(!is.null(file))
			dev.off()
	}
}


#-----------------------------------------------------------------------
#    Print function for an object of class bcrm
# -----------------------------------
print.bcrm<-function(x,...){
	cat(" Estimation method: ",x$method,"\n")
	ff.txt<-switch(x$ff
		,ht="Hyperbolic Tangent"
		,logit1="1-parameter logistic"
		,power="1-parameter power"
		,logit2="Two-parameter logistic")

	cat("\n Model: ",ff.txt,"\n")

	pa.txt<-switch(x$prior.alpha[[1]]
		,"1"=paste("Gamma( Shape:",x$prior.alpha[[2]],", Scale:",x$prior.alpha[[3]],")",sep="")
		,"2"=paste("Uniform(",x$prior.alpha[[2]],", ",x$prior.alpha[[3]],")",sep="")
		,"3"=paste("Lognormal( Mean:",x$prior.alpha[[2]],", Variance:",x$prior.alpha[[3]],")",sep="")
		,"4"=paste("Log Multivariate Normal"))
	cat("\n Prior: ",pa.txt,"\n")
	if(x$prior.alpha[[1]]==4){
		cat("Mean Vector: \n")
		print(x$prior.alpha[[2]])
		cat("\nVariance-Covariance Matrix: \n")
		print(x$prior.alpha[[3]])
	}
	tab1<-x$sdose
	names(tab1)<-x$dose
	cat("\n Standardised doses (skeleton): \n")
	print(tab1)

	if(x$constrain) { 
		if(is.null(x$dose)){
			cat("\n Modified (constrained) CRM used, starting dose level: ",x$start,"\n") 
		} else {
			cat("\n Modified (constrained) CRM used, starting dose: ",x$dose[x$start],"\n") 
		}
	} else { cat("\n Unmodified (unconstrained) CRM used \n") }
	
	if(!is.null(x$loss)){
		cat("\n Loss function given intervals of toxicity used to select next dose.")
		tab.lf<-x$loss
		names(tab.lf)<-levels(cut(0,breaks=c(0,x$tox.cutpoints,1)))
		cat("\n Loss function: \n")
		print(tab.lf)
	} else if(x$pointest=="plugin"){
		cat("\n Plug-in estimate of probability of toxicity used to select next dose \n")
	} else if(x$pointest=="mean"){
		cat("\n Posterior mean estimate of probability of toxicity used to select next dose \n")
	} else { 
		cat("\n",100*x$pointest,"percentile of (standardised) MTD distribution used to select next dose")
		cat("\n",100*x$pointest,"percentile is:",x$ndose[[length(x$ndose)]]$target,"\n")
	}	

	tab<-rbind(x$tox+x$notox,x$tox)
	rownames(tab)<-c("n","Toxicities")
	colnames(tab)<-x$dose
	names(dimnames(tab))<-c("","Doses")
	cat("\n Toxicities observed: \n")
	print(tab)

	if(x$method %in% c("exact","exact.sim") & x$ff=="logit2"){
		tab2a<-signif(rbind(x$ndose[[length(x$ndose)]]$est),3)
		rownames(tab2a)<-c("Point estimate")
	} else {
		tab2a<-signif(rbind(x$ndose[[length(x$ndose)]]$mean,x$ndose[[length(x$ndose)]]$sd,x$ndose[[length(x$ndose)]]$quantiles["50%",]),3)
		rownames(tab2a)<-c("Mean","SD","Median")
		tab2b<-signif(x$ndose[[length(x$ndose)]]$quantiles,3)
		colnames(tab2b)<-x$dose
		names(dimnames(tab2b))<-c("Quantiles","Doses")
	}
	colnames(tab2a)<-x$dose
	names(dimnames(tab2a))<-c("","Doses")
	cat("\n Posterior estimates of toxicity: \n")
	print(tab2a)
	if(!(x$method %in% c("exact","exact.sim") & x$ff=="logit2")){
		print(tab2b)
	}
	if(!is.null(x$loss)){
		tab3<-rbind(x$ndose[[length(x$ndose)]]$est)
		colnames(tab3)<-x$dose
		cat("\n Posterior expected loss at each dose: \n")		
		print(tab3)
	} else if(x$pointest=="plugin"){
		tab3<-signif(rbind(x$ndose[[length(x$ndose)]]$est),3)
		colnames(tab3)<-x$dose
		cat("\n Plug-in estimates of toxicity: \n")
		print(tab3)
	}
	if(is.null(x$dose)){
		cat("\n Next recommended dose level: ",x$ndose[[length(x$ndose)]]$ndose,"\n")
	} else {
		cat("\n Next recommended dose: ",x$dose[x$ndose[[length(x$ndose)]]$ndose],"\n")
	}
}

#-----------------------------------------------------------------------
#    Print function for an object of class bcrm.sim
# -----------------------------------
print.bcrm.sim<-function(x,tox.cutpoints=NULL,trajectories=FALSE,...){
	if(trajectories){
		## sample size
		n<-sapply(x,function(i){length(i$alldoses)})
		dose.mat<-outcome.mat<-matrix(NA,nrow=length(x),ncol=max(n))
		for(i in 1:length(x)){
			dose.mat[i,1:n[i]]<-x[[i]]$alldoses
			outcome.mat[i,1:n[i]]<-x[[i]]$alltox
		}
		return(list(doses=dose.mat,outcomes=outcome.mat))
	} else {
		# average sample size
		n.average<-mean(sapply(x,function(i){length(i$alldoses)}))
		n.min<-min(sapply(x,function(i){length(i$alldoses)}))
		n.max<-max(sapply(x,function(i){length(i$alldoses)}))
		if(n.min==n.max){
			tab0<-cbind(n.average)
			rownames(tab0)<-"Sample size"
			colnames(tab0)<-""
		} else {
			tab0<-cbind(n.average,n.min,n.max)
			rownames(tab0)<-"Sample size"
			colnames(tab0)<-c("Mean","Minimum","Maximum")
		}
		exp<-apply(sapply(x,function(i){(i$tox+i$notox)/sum(i$tox+i$notox)}),1,mean)
		rec<-prop.table(table(factor(sapply(x,function(i){i$ndose[[length(i$ndose)]]$ndose}),levels=1:length(x[[1]]$tox))))
		tab<-signif(rbind(exp,rec),3)
		rownames(tab)<-c("Experimentation proportion","Recommendation proportion")
		colnames(tab)<-x[[1]]$dose
		names(dimnames(tab))<-c("","Doses")
		if(is.null(tox.cutpoints)){
			tox.cutpoints<-seq(0,1,by=0.2)
		} else {
			tox.cutpoints<-unique(c(0,tox.cutpoints,1))
		}
		exp.tox<-prop.table(table(cut(unlist(sapply(x,function(i){rep(i$truep,(i$tox+i$notox))},simplify=FALSE)),tox.cutpoints,include.lowest=T)))
		rec.tox<-prop.table(table(cut(sapply(x,function(i){i$truep[i$ndose[[length(i$ndose)]]$ndose]}),tox.cutpoints,include.lowest=T)))
		tab2<-signif(rbind(exp.tox,rec.tox),3)
		rownames(tab2)<-c("Experimentation proportion","Recommendation proportion")
		names(dimnames(tab2))<-c("","Probability of DLT")
		cat("Operating characteristics based on ",length(x)," simulations: \n \n")
		print(tab0)
		cat("\n")
		print(tab)
		cat("\n")
		print(tab2)
	}
}

#-----------------------------------------------------------------------
## HT Gamma
## WinBUGS code for a hyperbolic tangent model with Gamma prior on alpha 
#-----------------------------------------------------------------------
HTGamma<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))), alpha)
	}
	alpha ~ dgamma(p1, p2)
}

#-----------------------------------------------------------------------
## HT Unif
## WinBUGS code for a hyperbolic tangent model with Uniform prior on alpha 
#-----------------------------------------------------------------------
HTUnif<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))), alpha)
	}
	alpha ~ dunif(p1, p2)
}
#-----------------------------------------------------------------------
## HT Lognormal
## WinBUGS code for a hyperbolic tangent model with Lognormal prior on alpha 
#-----------------------------------------------------------------------
HTLogNormal<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))), alpha)
	}
	prec<-1/p2
	alpha ~ dlnorm(p1,prec)
}

#-----------------------------------------------------------------------
## Logistic Gamma
## WinBUGS code for a one-parameter logistic model with Gamma prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticGamma<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
	}
	alpha ~ dgamma(p1, p2)
}

#-----------------------------------------------------------------------
## Logistic Uniform
## WinBUGS code for a one-parameter logistic model with Uniform prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticUnif<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
	}
	alpha ~ dunif(p1, p2)
}

#-----------------------------------------------------------------------
## Logistic Lognormal
## WinBUGS code for a one-parameter logistic model with Lognormal prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticLogNormal<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
	}
	prec<-1/p2
	alpha ~ dlnorm(p1, prec)
}

#-----------------------------------------------------------------------
## Power Gamma
## WinBUGS code for a power model with Gamma prior on alpha
#-----------------------------------------------------------------------
PowerGamma<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow(d[i], alpha)
	}
	alpha ~ dgamma(p1, p2)
}

#-----------------------------------------------------------------------
## Power Uniform
## WinBUGS code for a power model with Uniform prior on alpha
#-----------------------------------------------------------------------
PowerUnif<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow(d[i], alpha)
	}
	alpha ~ dunif(p1, p2)
}
## Power LogNormal
## WinBUGS code for a power model with LogNormal prior on alpha
#-----------------------------------------------------------------------
PowerLogNormal<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		p[i] <- pow(d[i], alpha)
	}
	prec<-1/p2
	alpha ~ dlnorm(p1, prec)
}
#-----------------------------------------------------------------------
## 2-parameter Logistic Bivariate Lognormal
## WinBUGS code for a two-parameter logistic model with Bivariate Lognormal prior on alpha[1] (intercept) and alpha[2] (slope of logistic function) 
#-----------------------------------------------------------------------
TwoPLogisticLogNormal<-function(){
	for (i in 1:N1){
		s[i] ~ dbin(p[i], n[i])		
		logit(p[i]) <- log.alpha[1] + alpha[2] * d[i]
	}
	alpha[1]<-exp(log.alpha[1])
	alpha[2]<-exp(log.alpha[2])
	Omega[1:2,1:2]<-inverse(p2[,])
	log.alpha[1:2] ~ dmnorm(p1[], Omega[,])
}

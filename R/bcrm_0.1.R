## Bayesian CRM - extending original code of J. Jack Lee and Nan Chen, Department of Biostatistics, the University of Texas M. D. Anderson Cancer Center
## Now using exact inference, R2WinBUGS and BRugs since BRugs is currently unavailable on CRAN 

### MJS 01/09/11

# ----------------------------------------------------------------------
# 	bcrm. Conduct a Bayesian CRM trial simulating outcomes from a true model
# 	N	     --> size of trial
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#	p.tox0     --> prior probabilities of response at dose levels
#     sdose      --> standardised dose levels (given if p.tox0 is missing)
#    dose		--> optional vector of actual doses (for plotting)
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
# 	pointest   --> Which summary estimate of the posterior distribution should be used to choose next dose, "plugin" (default), "mean" or a numerical quantile between 0 and 1 (e.g. 0.5)
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
# ----------------------------------------------------------------------
bcrm<-function(N,tox=NULL,notox=NULL,p.tox0=NULL,sdose=NULL,dose=NULL,ff,prior.alpha,cohort=3,target.tox,constrain=TRUE,sdose.calculate="mean",pointest="plugin",start=NULL,simulate=FALSE,nsims=1,truep=NULL,
	method="exact",burnin.itr=2000,production.itr=2000,bugs.directory="c:/Program Files/WinBUGS14/",plot=FALSE,file=NULL){

	if(missing(p.tox0) & missing(sdose)) stop("Either p.tox0 or sdose must be specified")
	if(!missing(p.tox0) & !missing(sdose)) stop("Only one of p.tox0 and sdose must be specified")
	if(sdose.calculate!="mean" & sdose.calculate!="median") stop("sdose.calculate must be either `mean' or `median'")
	if((is.character(pointest) & pointest!="mean" & pointest!="plugin") | is.numeric(pointest) & (pointest<0 | pointest>1)) stop("pointest must be either `plugin', `mean' or an EWOC feasibility quantile between 0 and 1")
	if(is.numeric(pointest) & method=="exact") stop("EWOC design must be fitted using MCMC methods")
	if(simulate & is.null(truep)) stop("truep must be specified if simulating data")
	if(simulate){ 
		plot<-FALSE
		results<-list()
	}
	if(!(method %in% c("exact","BRugs","R2WinBUGS"))) stop("method must be either `exact', `BRugs' or `R2WinBUGS'")
	if(constrain & is.null(start)) stop("A starting/current dose level must be specified using `start' if constrain==TRUE")
	## Changing functional form to numeric 
  ff<-ifelse(ff=="ht",1,ifelse(ff=="logit1",2,ifelse(ff=="power",3,ifelse(ff=="logit2",4,stop("ff must be one of `ht', `logit1', `power' or `logit2'")))))
  if(ff==4 & (length(prior.alpha[[2]])<2 | length(prior.alpha[[3]])<2)) stop("second and third components of `prior.alpha' must be vectors of size 2")

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

	k<-length(sdose)
	sim<-1	
	
	while(sim<=nsims){
		## Prior
		if(is.null(tox) & is.null(notox)){
			new.tox<- rep(0,k)
			new.notox<-rep(0,k)
			current<-start-1
			alpha<-switch(method
			,BRugs=getprior(prior.alpha, 10000)
			,R2WinBUGS=getprior(prior.alpha,10000)
			,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha))
		} else { 
			new.tox<-tox
			new.notox<-notox
			current<-start 
			alpha<-switch(method
			,BRugs=Posterior.BRugs(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
			,R2WinBUGS=Posterior.R2WinBUGS(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr,bugs.directory)
			,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
			)
		}
		ncurrent<-sum(new.tox+new.notox)
		alldoses<-NULL
		ndose<-if(method!="exact") nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current)
				else nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
		if(plot){
			dose.response.plot(new.tox,new.notox,dose,sdose,ff,target.tox,ndose)
			if(!is.null(file))
			ggsave(paste(file,ncurrent,".pdf",sep=""))
		}
		while(ncurrent<N){
			current<-ndose[[1]]
			alldoses<-c(alldoses,rep(current,cohort))
			ncurrent<-ncurrent+cohort
			if(!simulate){
				interact<-crm.interactive(new.tox,new.notox,ncurrent,cohort,ndose)
				if(interact$bk==TRUE){
					results<-list(dose=dose,sdose=sdose,tox=new.tox,notox=new.notox,alpha=alpha,ndose=ndose,constrain=constrain,start=start,target.tox=target.tox,ff=ff,method=method,pointest=pointest,prior.alpha=prior.alpha,alldoses=alldoses)
					class(results)<-"bcrm"
					return(results)
				}
				new.tox<-interact$tox
				new.notox<-interact$notox
			} else {	y<-rbinom(1,cohort,truep[ndose[[1]]]) 
				new.notox[ndose[[1]]]<-new.notox[ndose[[1]]]+(cohort-y)
				new.tox[ndose[[1]]]<-new.tox[ndose[[1]]]+y
			}
			alpha<-switch(method
				,BRugs=Posterior.BRugs(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
				,R2WinBUGS=Posterior.R2WinBUGS(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr,bugs.directory)
				,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha))
			ndose<-if(method!="exact") nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current)
				else nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
			if(plot){
				dose.response.plot(new.tox,new.notox,dose,sdose,ff,target.tox,ndose)
				if(!is.null(file))
				ggsave(paste(file,ncurrent,".pdf",sep=""))
			}
		}
		if(!simulate){
			results<-list(dose=dose,sdose=sdose,tox=new.tox,notox=new.notox,alpha=alpha,ndose=ndose,constrain=constrain,start=start,target.tox=target.tox,ff=ff,method=method,pointest=pointest,prior.alpha=prior.alpha,alldoses=alldoses)
			class(results)<-"bcrm"
	    	} else {
			results[[sim]]<-list(dose=dose,sdose=sdose,tox=new.tox,notox=new.notox,ndose=ndose,constrain=constrain,target.tox=target.tox,start=start,method=method,pointest=pointest,prior.alpha=prior.alpha,alldoses=alldoses)
			class(results)<-"bcrm.sim"
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
# ----------------------------------------------------------------------
crm.interactive<-function(tox,notox,ncurrent,cohort,ndose){
     k <- length(tox)
		repeat {
     		  cat("\n\n RECOMMENDED DOSE LEVEL FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent, "IS:", ndose[[1]])
	            ans <- get.dose(k,ncurrent,cohort)     
     	       if (ans==-2)
          	      ans <- ndose[[1]]            
	            if (ans==0)
              	  cat("\n\n EXIT AND RETURN THE RESULTS SO FAR?")            
	            if (ans<1) {
     	           cat("\n\n DO YOU REALLY WANT TO EXIT ? (Y/N)  ")
          	      yn <- readline()
              	 	 if (yn=="y" || yn=="Y") break
              		  else                    next
	     		 }  
			# y.j is the number of toxicities from the treatment of
		     # cohort patients at dose level ndose
     		cat("\n ENTER NUMBER OF TOXICITIES FOR PATIENTS ",ncurrent-cohort+1," to ",ncurrent,":\n")
		     y <- get.answer(cohort)               
     		# give the user a last chance to modify the treatment and outcome
		     cat("\n\n\t\t ENTERED VALUES:")
     		cat("\n DOSE LEVEL ...", ans)
	     	cat("\n NO. OF TOXICIES ....",y)
	     	cat("\n PRESS `RETURN' IF OK, OR ANY OTHER KEY TO ENTER NEW VALUES ")
		     key <- readline()
     		if (nchar(key)==0) break
		}
          if (ans==0)  return(list(bk=TRUE))
		notox[ans]<-notox[ans]+(cohort-y)
		tox[ans]<-tox[ans]+y
		return(list(tox=tox,notox=notox,bk=FALSE))
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
get.dose <- function( n ,ncurrent,cohort) {
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
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter logistic
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
        if (ff == 1)
            HTGamma
        else if (ff == 2)
            LogisticGamma
        else if (ff == 3)
            PowerGamma
	   else stop("Functional form not currently available with specified prior distribution")
    }
    else if (prior.alpha[[1]] == 2)
    {
        if (ff == 1)
            HTUnif
        else if (ff == 2)
            LogisticUnif
        else if (ff == 3)
            PowerUnif
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 3)
    {
        if (ff == 1)
            HTLogNormal
        else if (ff == 2)
            LogisticLogNormal
        else if (ff == 3)
            PowerLogNormal
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 4)
    {
        if (ff == 4)
		TwoPLogisticLogNormal            
	   else stop("Functional form not currently available with specified prior distribution")
    }    
	path.model<-file.path(tempdir(), "model.file.txt")
	path.data<-file.path(tempdir(), "data.file.txt")
	write.model(model.file,path.model)
	bugsData(mydata,path.data)
	BRugsFit(path.model,path.data, inits=NULL, numChains = 2, parametersToSave="alpha",
	    	nBurnin = burnin.itr, nIter = production.itr/2, BRugsVerbose = FALSE,DIC=F)
	if(ff==4){
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
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter logistic
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
        if (ff == 1)
            HTGamma
        else if (ff == 2)
            LogisticGamma
        else if (ff == 3)
            PowerGamma
	   else stop("Functional form not currently available with specified prior distribution")
    }
    else if (prior.alpha[[1]] == 2)
    {
        if (ff == 1)
            HTUnif
        else if (ff == 2)
            LogisticUnif
        else if (ff == 3)
            PowerUnif
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 3)
    {
        if (ff == 1)
            HTLogNormal
        else if (ff == 2)
            LogisticLogNormal
        else if (ff == 3)
            PowerLogNormal
	   else stop("Functional form not currently available with specified prior distribution")
    }        
    else if (prior.alpha[[1]] == 4)
    {
        if (ff == 4)
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
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter logistic
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
			,"4"=function(alpha,prior.alpha){1/(alpha[1]*alpha[2])*dmvnorm(log(alpha),mean=prior.alpha[[2]],sigma=prior.alpha[[3]])})

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
		alpha.quantiles<-sapply(c(0.025,0.25,0.5,0.75,0.975),function(q){optimize(fn,c(0,10*alpha.mean),q=q, tol = .Machine$double.eps^0.50)$minimum})
		dose.quantiles<-if(ff==2){ sapply(sdose,function(d){wmodel(d,sort(alpha.quantiles,TRUE))})
		} else {sapply(sdose,function(d){wmodel(d,sort(alpha.quantiles,TRUE))})} 
		rownames(dose.quantiles)<-c("2.5%","25%","50%","75%","97.5%")
	} else {
		stop("Exact method not yet implemented for 2-parameter model")
	}
	return(list(alpha.mean=alpha.mean,dose.mean=dose.mean,dose.sd=dose.sd,dose.quantiles=dose.quantiles))
}


# ----------------------------------------------------------------------
#     find the new dose and posteria mean of alpha
#
#     samples.alpha  --> sampled data of variable(s) alpha
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#	constrain	  --> TRUE or FALSE
#	pointest	  --> "plugin", "mean" or a quantile
#	current	  --> current dose level (of last treated cohort) - 0 if this is first cohort
#
# ----------------------------------------------------------------------
nextdose<-function(samples.alpha,sdose,ff,target.tox,constrain,pointest,current){
	f<-which.f(ff)
	k<-length(sdose)
	if(pointest=="plugin"){
		mean.alpha<-apply(as.matrix(samples.alpha),2,mean)
		if(length(mean.alpha)>1) mean.alpha<-t(as.matrix(mean.alpha))
	} 
	samples.sdose<-sapply(sdose,function(x){f(x,samples.alpha)})
	if(is.numeric(pointest)){
		mtd<-find.x(ff,ptox=target.tox,alpha=samples.alpha)
		target<-quantile(mtd,pointest)
	} else {
		target<-target.tox
	}
	est<-if(pointest=="plugin"){ f(sdose,mean.alpha)
	} else {if(pointest=="mean"){apply(samples.sdose,2,mean)
		} else { sdose }}
	mean<-apply(samples.sdose,2,mean)
	median<-apply(samples.sdose,2,median)
	sd<-apply(samples.sdose,2,sd)
	quantiles<-apply(samples.sdose,2,quantile,c(0.025,0.25,0.50,0.75,0.975))
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
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter Logistic
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


# -----------------------------------------------------------------------
#     Return the calculation models
#     ff       --> integer value:
#                    1 - hyperbolic tangent
#                    2 - logistic
#                    3 - Power
#				  4 - Two-parameter logistic
# -----------------------------------------------------------------------
which.f <- function( ff ) {
    return( switch(ff,
                    function(dose,alpha) {((tanh(dose)+1)/2)^alpha},
                    function(dose,alpha) {1/(1+exp(-3-alpha*dose))},
                    function(dose,alpha) {dose^alpha},
                    function(dose,alpha) {1/(1+exp(-log(alpha[,1])-alpha[,2]*dose))} ))
}

# ----------------------------------------------------------------------
#     Find the dose corresponding to a certain p(DLT)
#
#     ff         --> functional form for the probability density
#     ptox       --> desired probability of DLT
#     alpha      --> the parameter(s) of the model
# ----------------------------------------------------------------------
find.x <- function( ff, ptox, alpha ) {
    if ( ff == 1 ) {
	  x<-atanh(2*ptox^(1/alpha)-1)
    }
    else if ( ff == 2 )
        x <- (qlogis(ptox)-3)/alpha
    else if ( ff == 3)
        x <- exp(log(ptox)/alpha)
	else if (ff==4){
	   if(is.vector(alpha))	alpha<-matrix(alpha,ncol=2)
	   x <- (qlogis(ptox)-log(alpha[,1]))/alpha[,2]
	}
	return( x )
}

#-----------------------------------------------------------------------
#     Dose-response plot 
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
# 	samples.alpha --> posterior sample of parameter(s) alpha
#	 dose	    --> vector of doses 
#     sdose       --> vector of standardised doses (used only if dose is not supplied)
#	ff		--> Functional form of the dose-response curve
#     target.tox  --> desired level of probability for toxicity
#	ndose	   --> summary information about next dose, together with posterior mean, sd and quantiles obtained from nextdose function
# -----------------------------------
dose.response.plot<-function(tox,notox,dose,sdose,ff,target.tox,ndose){
	dose<-if(is.null(dose)) sdose else dose
	f<-which.f(ff)

	df<-data.frame(dose=dose,target.tox=target.tox,est=ndose$est,mean=ndose$mean,q2.5=ndose$quantiles["2.5%",],q25=ndose$quantiles["25%",],q50=ndose$quantiles["50%",],q75=ndose$quantiles["75%",],q97.5=ndose$quantiles["97.5%",])
	df2<-data.frame(dose=factor(c(rep(dose,tox),rep(dose,notox)),levels=dose),Outcome=c(rep("DLT",sum(tox)),rep("No DLT",sum(notox))))

	a<-ggplot()+geom_errorbar(aes(x=dose,ymin=q2.5,ymax=q97.5),colour="red",data=df)+geom_crossbar(aes(x=dose,y=q50,ymin=q25,ymax=q75),data=df,fill="red")+
		geom_hline(aes(yintercept=target.tox),data=df,col=4,linetype=2)+geom_point(aes(x=dose,y=est),data=df[ndose[[1]],],size=4,col=4,shape=9)+xlab("Dose")+ylab("Probability of DLT")+ylim(0,1)+
		opts(title="Posterior p(DLT) quantiles: 2.5%, 25%, 50%, 75%, 97.5%")
	b<-if(nrow(df2)!=0) {ggplot()+geom_bar(aes(x=dose,fill=Outcome),data=df2)+xlab("Dose")+ylab("Number")
		} else { NULL	}

	grid.newpage()
	pushViewport(viewport(layout=grid.layout(2,1)))
	vplayout<-function(x,y)	viewport(layout.pos.row=x,layout.pos.col=y)
	print(a,vp=vplayout(1,1))	
	print(b,vp=vplayout(2,1))	
}

#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm
# -----------------------------------
plot.bcrm<-function(x,...){
	dose.response.plot(x$tox,x$notox,x$dose,x$sdose,x$ff,x$target.tox,x$ndose)
}



#-----------------------------------------------------------------------
#    Print function for an object of class bcrm
# -----------------------------------
print.bcrm<-function(x,...){
	cat(" Estimation method: ",x$method,"\n")
	ff.txt<-switch(x$ff
		,"1"="Hyperbolic Tangent"
		,"2"="1-parameter logistic"
		,"3"="1-parameter power"
		,"4"="Two-parameter logistic")

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

	if(x$constrain) { cat("\n Modified (constrained) CRM used, starting dose: ",x$start,"\n") 
		} else { cat("\n Unmodified (unconstrained) CRM used \n") }
	
	if(x$pointest=="plugin"){
		cat("\n Plug-in estimate of probability of toxicity used to select next dose \n")
	} else if(x$pointest=="mean"){
		cat("\n Posterior mean estimate of probability of toxicity used to select next dose \n")
	} else { 
		cat("\n",100*x$pointest,"percentile of (standardised) MTD distribution used to select next dose")
		cat("\n",100*x$pointest,"percentile is:",x$ndose$target,"\n")
	}	

	tab<-rbind(x$tox+x$notox,x$tox)
	rownames(tab)<-c("n","Toxicities")
	colnames(tab)<-x$dose
	names(dimnames(tab))<-c("","Doses")
	cat("\n Toxicities observed: \n")
	print(tab)
	
	tab2a<-rbind(x$ndose$mean,x$ndose$sd,x$ndose$quantiles["50%",])
	rownames(tab2a)<-c("Mean","SD","Median")
	colnames(tab2a)<-x$dose
	names(dimnames(tab2a))<-c("","Doses")
	tab2b<-x$ndose$quantiles
	colnames(tab2b)<-x$dose
	names(dimnames(tab2b))<-c("Quantiles","Doses")
	cat("\n Posterior estimates of toxicity: \n")
	print(tab2a)
	print(tab2b)
	if(x$pointest=="plugin"){
		tab3<-rbind(x$ndose$est)
		colnames(tab3)<-x$dose
		cat("\n Plug-in estimates of toxicity: \n")
		print(tab3)
	}
	cat("\n Next recommended dose level: ",x$ndose$ndose," (",x$dose[x$ndose$ndose],")\n")
	
	
	
}

#-----------------------------------------------------------------------
#    Print function for an object of class bcrm.sim
# -----------------------------------
print.bcrm.sim<-function(x,...){
	n<-apply(sapply(x,function(i){(i$tox+i$notox)/sum(i$tox+i$notox)}),1,mean)
	rec<-prop.table(table(factor(sapply(x,function(i){i$ndose}),levels=1:length(x[[1]]$tox))))
	tab<-rbind(n,rec)
	rownames(tab)<-c("Experimentation proportion","Recommendation proportion")
	colnames(tab)<-x[[1]]$dose
	names(dimnames(tab))<-c("","Doses")
	cat("Operating characteristics based on ",length(x)," simulations: \n")
	print(tab)
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

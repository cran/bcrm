simFun <- function(X, stop, ndose, sdose, dose, constrain, start, ff,
                   cohort, method, pointest, tox.cutpoints, loss,
                   burn, iter, quantiles, target.tox,
                   prior.alpha, truep, only.below, quietly=10, bugs.directory){

  k <- length(sdose)
  
  new.tox <- new.notox <- rep(0, k)
  current <- start - 1
  first<-FALSE
  ncurrent <- sum(new.tox + new.notox)
  newdata <- data.frame(patient=NULL, dose=NULL, tox=NULL)
  
  stopped <- stop.check(stop, target.tox, ncurrent, ndose, new.tox, new.notox,
                        simulate = TRUE)
  ndose <- stopped$ndose ## update ndose in case no doses are deemed safe
  
  results <- list(dose=dose, sdose=sdose, tox=new.tox,
                  notox=new.notox, ndose=list(ndose),
                  constrain=constrain, start=start,
                  target.tox=target.tox, ff=ff,
                  method=method, pointest=pointest,
                  tox.cutpoints=tox.cutpoints,
                  loss=loss, prior.alpha=prior.alpha,
                  truep=truep, only.below=only.below)
  
  while(!stopped$stop){
    current <- ndose[[1]]
    ncurrent <- ncurrent+cohort
    
    y <- rbinom(cohort, 1, truep[ndose[[1]]]) 
    new.notox[ndose[[1]]] <- new.notox[ndose[[1]]] + (cohort - sum(y))
    new.tox[ndose[[1]]] <- new.tox[ndose[[1]]] + sum(y)
    
    currentdata <- data.frame(patient = (ncurrent - cohort + 1):ncurrent,
                              dose = rep(current, cohort),
                              tox = y)
    newdata <- rbind(newdata, currentdata)
    
    alpha <- switch(method
                    , rjags=Posterior.rjags(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burn,  iter)
                    , BRugs=Posterior.BRugs(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burn,  iter)
                    , R2WinBUGS=Posterior.R2WinBUGS(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burn, iter, bugs.directory)
                    , exact=Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
                    , exact.sim=Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
    )
    
    ndose <- switch(method
                    , rjags=nextdose(alpha, sdose, ff, target.tox, constrain, first, pointest, current, tox.cutpoints, loss, quantiles, only.below)
                    , BRugs=nextdose(alpha, sdose, ff, target.tox, constrain, first, pointest, current, tox.cutpoints, loss, quantiles, only.below)
                    , R2WinBUGS=nextdose(alpha, sdose, ff, target.tox, constrain, first, pointest, current, tox.cutpoints, loss, quantiles, only.below)
                    , exact=nextdose.exact(alpha, sdose, ff, target.tox, constrain, first, pointest, current, only.below)
                    , exact.sim=nextdose.exact.sim(alpha, sdose, ff, target.tox, constrain, first, pointest, current, only.below)
    )
    
    stopped <- stop.check(stop, target.tox, ncurrent, ndose, new.tox, new.notox,
                          simulate=TRUE)
    ndose <- stopped$ndose ## update ndose in case no doses are deemed safe
    
    results$tox <- new.tox
    results$notox <- new.notox
    results$ndose[[length(results$ndose)+1]] <- ndose
    results$data <- newdata
  } # Close while
  ## Once stopped run following code
  if(X %% quietly == 0){
    message("Simulated trial: ", X)
  }

  results
} # Close simFun
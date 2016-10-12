matrixordpls.boot <- function(data, model, ..., R = 500,
                           signChange = NULL,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("boot.ncpus", 1L),
                           dropInadmissible = FALSE,
                           stopOnError = FALSE,
                           extraFun = NULL){
  
  if(! requireNamespace("boot")) stop("matrixpls.boot requires the boot package")
  
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  
  
  model <- matrixpls:::parseModelToNativeFormat(model)
  
  data <- data[,rownames(model$reflective)]
  data <- as.matrix(data)
  
  arguments <- list(model = model, ...)
  
  # Prepare sign change corrections
  
  if(! is.null(signChange)){
    # Worig <- attr(matrixpls(stats::cov(data),...), "W")
    Worig <- attr(matrixpls::matrixpls(polycor::hetcor(data,std.err=FALSE,ML=FALSE)$correlations,...), "W")
    
    # Get the original weight function
    weightFun <- arguments[["weightFunction"]]
    if(is.null(weightFun)) weightFunction <- matrixpls:::weightFun.pls
    
    # Wrap inside sign change correction
    arguments[["weightFun"]] <- function(S, ...){
      Wrep <- weightFun(S, ...)
      W <- matrixpls:::signChange(Worig, Wrep)
      W <- matrixpls:::scaleWeights(S, W)
      attributes(W) <- attributes(Wrep)
      W
    }
    
  }
  
  # Bootstrap
  
  boot.out <- boot::boot(data,
                         function(data, indices, ...){
                           
                           # S <- stats::cov(data[indices,])
                           S=polycor::hetcor(data[indices,],std.err=FALSE,ML=FALSE)$correlations #polychoric correlation is used instead of the BP correlation
                           arguments <- c(list(S),arguments)
                           
                           if(stopOnError){
                             boot.rep <- do.call(marixpls::matrixpls, arguments)
                           }
                           else{
                             tryCatch(
                               boot.rep <- do.call(matrixpls::matrixpls, arguments)
                             )
                           }
                           
                           
                           # Add additional statistics
                           
                           if(!is.null(extraFun)){
                             a <- attributes(boot.rep)
                             boot.rep <-c(boot.rep,extraFun(boot.rep))
                             a$names <- names(boot.rep)
                             attributes(boot.rep) <- a
                           }
                           
                           # Add convergence status as the last statistic. This will be removed
                           # later
                           
                           boot.rep[length(boot.rep)+1] <- matrixpls:::convergenceStatus(boot.rep)
                           
                           # If the indices are not sorted, then this is not the original sample
                           # and we can safely omit all attributes to save memory. 
                           
                           if(is.unsorted(indices)){
                             attributes(boot.rep) <- NULL
                           }
                           
                           boot.rep
                           
                         },
                         R, parallel = parallel, ncpus = ncpus, ...)
  
  
  # Store the convergence status frequencies
  boot.out$convergence <- table(boot.out$t[,1])
  
  # Clean inadmisibles
  
  i <- length(boot.out$t0)
  
  if(dropInadmissible){
    boot.out$t <- boot.out$t[boot.out$t[,i]==0,]
    boot.out$R <- nrow(boot.out$t)
  }
  
  # Remove the convergence status from the estimates
  a <- attributes(boot.out$t0)
  a$names <- a$names[-length(boot.out$t0)]
  boot.out$t <- boot.out$t[,-i]
  boot.out$t0 <- boot.out$t0[-i]
  attributes(boot.out$t0) <- a
  
  # Set class and return
  class(boot.out) <- c("matrixplsboot", class(boot.out))
  boot.out
}

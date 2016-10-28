source("likelihoods.R") 
generateTestData <- function(forms, individual_params, shared_params, n=1){
   means <- sapply(forms, eval, c(as.list(individual_params), as.list(shared_params)))
   counts <- rpois(n=length(means)*n, lambda=means)
   data.frame(condition=rep(names(forms), n), count=counts)
}

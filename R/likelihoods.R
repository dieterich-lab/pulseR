#from pryr package
substitute_q <- function(x, env)
{
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

# get matrix for samples
sample_means <- function(evaled_forms, norm_factors){
  evaled_forms %*% norm_factors
}

# universal likehood
ll <- function(par, namesToOptimise, pd, byOne=FALSE) {
  pattern <- par[namesToOptimise]
  if (byOne)
    pattern <- lapply(pattern, `[[`, 1)
  par[namesToOptimise] <- NULL
  evalCall <- as.call(c(cbind, pd$formulas))
  norms <- getNorms(pd, par$normFactors)
  if (!byOne)
    getPars <- quote(relist(x, pattern))
  else 
    getPars <- quote(as.list(x))
  f <- bquote(
    function(x, counts, fixedPars = par) {
      fixedPars[namesToOptimise] <- .(getPars)
      evaledForms <- eval(evalCall, fixedPars)
      means <- sample_means(evaledForms, norms)
      -sum(stats::dnbinom(
        counts, mu = means, size = fixedPars$size, log = TRUE))
    }, list(getPars = getPars)) 
  eval(f)
}

getNorms <- function(pd, normFactors = NULL) {
  m <- matrix(0,
           ncol = length(pd$formulaIndexes),
           nrow = length(pd$formulas))
  indexes <- do.call(rbind,lapply(seq_along(pd$formulaIndexes),
         function(i){
           cbind(pd$formulaIndexes[[i]],i)
         }))
  norms <- unlist(pd$depthNormalisation)
  if (!is.null(normFactors)) {
    norms <- norms * unlist(normFactors)[unlist(pd$interSampleIndexes)]
  }
  m[indexes] <- norms
  m
}

# likelihood for norm factors
llnormFactors <- function(par, pd) {
  evaledForms <- eval(as.call(c(cbind, pd$formulas)), par)
  function(x, counts) {
    norms <- getNorms(pd, c(1,x))
    means <- sample_means(evaledForms, norms)
    -sum(stats::dnbinom(counts, mu = means, size = par$size, log = TRUE))
  }
}

totalll <- function(par, pd) {
  function(x, counts) {
    x <- relist(x, par)
    evaledForms <- eval(as.call(c(cbind, pd$formulas)), par)
    norms <- getNorms(pd, c(1, x$normFactors))
    means <- sample_means(evaledForms,  norms)
    - sum(stats::dnbinom( counts,mu = means, size = x$size, log = TRUE))
  }
}

#' Calculate expected read numbers for the raw data 
#'
#' @param par estimated parameters from \link{fitModel}
#' @param pulseData a \link{PulseData} object.
#'
#' @return a list(preditions=,llog=), where
#'   predictions is a matrix of the same dimension 
#'   as the raw counts in pulseData$count_data;
#'   llog is a matrix with logarithms of likelihood for the given raw counts.
#' @export
#'
predictExpression <- function(pulseData, par) {
  evaledForms <- eval(as.call(c(cbind, pulseData$formulas)), par)
  norms <- getNorms(pulseData, par$normFactors)
  means <- sample_means(evaledForms, norms)
  llog <- stats::dnbinom(counts, mu = means, size = par$size, log = TRUE)
  list(predictions = means, llog = llog)
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}


evaluateLikelihood <- function(pulseData, par ) {
  objective <- ll_dispersion(pulseData, par)
  objective(par$size)
}

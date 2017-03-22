#from pryr package
substitute_q <- function(x, env)
{
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

#' Get matrix for samples
#' 
#' Means are defined as a linear combination of the formulas, and
#' we use matrix representation of this relations via norm_factors matrix.
#' 
#' For example, if there are pull-down and flow-through fractions, 
#' which have 10% contamination by each other, and their mean read numbers
#' are 10 and 20 respectevily, the computation may look like
#'         |.9 .1|
#' |10 20|x|     |
#'         |.1 .9|
#' for one pull-down and one flow-through samples.
#' 
#' @param evaled_forms a numeric vector of formulas evaluated for every
#' condition
#' @param norm_factors a matrix n x k, where k is the number of samples
#'
#' @return a vector of length equal to the sample number
#'
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

#' Calculates mean read number estimations 
#'
#' @param par estimated parameters from \link{fitModel}
#' @param pulseData a \link{PulseData} object.
#'
#' @return a named list:
#'   - preditions, a matrix of the same dimension as of the raw counts 
#'   - llog, a matrix with logarithms of likelihood for the given raw counts.
#' @export
#'
#' @examples 
#' \dontrun{
#' # Plot expected values vs the raw counts.
#' # Let res is the return of fitModel and pd is a PulseData object
#' pr <- predictExpression(pd, res)
#' plot(y = pr$predictions, x = pd$counts, xlab = "raw", ylab = "fitted")
#' }
predictExpression <- function(pulseData, par) {
  evaledForms <- eval(as.call(c(cbind, pulseData$formulas)), par)
  norms <- getNorms(pulseData, par$normFactors)
  means <- sample_means(evaledForms, norms)
  llog <- stats::dnbinom(
    pulseData$counts, mu = means, size = par$size, log = TRUE)
  list(predictions = means, llog = llog)
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}


evaluateLikelihood <- function(pulseData, par) {
  evaledForms <- eval(as.call(c(cbind, pulseData$formulas)), par)
  norms <- getNorms(pulseData, par$normFactors)
  means <- sample_means(evaledForms, norms)
  llog <- stats::dnbinom(
    pulseData$counts, mu = means, size = par$size, log = TRUE)
  sum(llog)
}

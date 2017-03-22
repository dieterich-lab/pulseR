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
#' By multiplying evaluated formulas vector by norm_factors matrix, 
#' we compute a weighted sum of the formulas. 
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

#' Creates a likelihood function for given parameter names
#'
#' @param par a named list with parameter values (including the size parameter
#'  for the negative binomial distribution, see \code{\link{dnbinom}})
#' @param namesToOptimise names of not fixed parameters
#' @param pd a \code{\link{PulseData}} object
#' @param byOne logical. If TRUE, the created function works on the level of
#' a single gene/isoform. In this case, every gene-specific parameter is
#' representd by a single scalar value. Otherwise, gene-specific parameters
#' are represented by vectors of the length equal to gene number.
#'
#' @return a function with the folowing arguments:
#'   - x, a numeric vector of variable parameters, which is used to 
#'     calculate the likelihood function. The order of parameters is as in 
#'     unlist(par)
#'   - counts, a numeric matrix, or, if byOne is TRUE, a numeric vector with
#'     read counts for every sample. If byOne is FALSE (default),
#'     number of rows must be equal the number of genes.
#'   - fixedPars is a list with the rest of parameters which are needed for
#'     likelihood calculation. By default, the ones provided in the `par`
#'     argument to `ll` function are used.
#'  The created function returns a logarithm of the likelihood value
#'  calculated on the basis of the negative binomial distribution for the
#'  provided counts and parameters.
#' @export
#'
ll <- function(par, namesToOptimise, pd, byOne=FALSE) {
  # we use relist-unlist idiom in the likelihood implementation
  # pattern is needed for recontruction of the correct parameter list
  # from the flatten vector of parameters, which is used by 
  # optimisation function
  pattern <- par[namesToOptimise]
  # leave only first parameter values if creating a gene-specific likelihood
  if (byOne)
    pattern <- lapply(pattern, `[[`, 1)
  par[namesToOptimise] <- NULL
  evalCall <- as.call(c(cbind, pd$formulas))
  norms <- getNorms(pd, par$normFactors)
  # for better performance, we don't relist in gene-specific likelihood
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

#' Constructs a matrix of normalisation coefficients
#'
#' @param pd \code{\link{PulseData}} object
#' @param normFactors a list of normalisation factors if inter-fraction
#' normalisation is used (i.e. no spike-ins are provided). 
#' If NULL (default), only sequencing depth normalisation is used
#' on the basis of spike-ins counts. The structure of `normFactors` must be
#' same as of `pd$interSampleCoeffs`.
#'
#' @return a matrix of normalisation coefficients to use during calculation of 
#' mean read number in samples. The row number equals number of 
#' different formulas used in estimation of means 
#' (i.e. the same as in `pd$rawFormulas`). The columns correspond to the samples
#' in the count matrix of the PulseData object `pd`.
#'
getNorms <- function(pd, normFactors = NULL) {
  m <- matrix(0,
           ncol = length(pd$formulaIndexes),
           nrow = length(pd$formulas))
  # create indexes for location of normalisation coefficients for every
  # sample in the normalisation matrix
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

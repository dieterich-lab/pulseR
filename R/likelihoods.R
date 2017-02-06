
#from pryr package
substitute_q <- function(x, env)
{
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

makeVector <- function(forms) {
  as.call(c(c,forms))
}

getNormFactors <- function(pulseData, par) {
  norm_factors <- pulseData$norm_factors
  if (!is.null(par$fraction_factors)) {
    norm_factors <-
      norm_factors * c(1, par$fraction_factors)[as.integer(pulseData$fraction)]
  }
  norm_factors
}

getMeans <- function(formulas, par) {
  params <- c(par$individual_params, par$shared_params, par$known)
  means <- lapply(formulas, function(x) {
    eval(x, envir = params)
  })
  means <- do.call(cbind, means) 
  means
}

#' Create a likelihood function for gene-specific parameters
#' 
#' The values of shared parameters, \code{size} from \code{\link{dnbinom}} and
#' normalisation factors are taken from \code{par}. 
#'
#' @param pulseData PulseData object
#' @param par list; must have fields \code{shared_params}, \code{size},
#'  \code{individual_params}, \code{fraction_factors}.
#'
#' @return a function(params, counts), which returns a  log likelihood
#' for a given vector of individual parameters, which are ordered as in 
#' \code{par$individual_params}. 
#' @importFrom stats dnbinom
#'
ll_gene <- function(pulseData, par) {
  mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
  formulas <- pulseData$formulas
  if (!is.null(par$shared_params))
    formulas <- lapply(formulas, substitute_q, par$shared_params)
  means_vector <-  makeVector(formulas)
  param_names <- names(par$individual_params)
  norm_factors <- getNormFactors(pulseData, par)
  function(params, counts, known=NULL) {
    names(params) <- param_names
    mus <- eval(means_vector, as.list(c(params, known)))
    if(any(mus<=0)) return(Inf)
    lambdas <-  mus[mean_indexes]
    - sum(dnbinom(
      x    = counts,
      mu   = lambdas * norm_factors,
      log  = TRUE,
      size = par$size
    ))
  }
}

#' Create a likelihood function for shared parameters
#' 
#' The values of gene-specific parameters, \code{size} from
#'  \code{\link{dnbinom}} and
#' normalisation factors are taken from \code{par}. 
#' @inheritParams ll_gene
#'
#' @return a function(params, counts), which returns a  log likelihood
#' for a given vector of shared parameters, which are ordered as in 
#' \code{par$shared_params}.
#' @importFrom stats dnbinom
#'
ll_shared_params <- function(pulseData, par) {
  shared_param_names <- names(par$shared_params)
  norm_factors <- getNormFactors(pulseData, par)
  function(shared_params) {
    names(shared_params) <- shared_param_names
    par$shared_params <- shared_params
    means <- getMeans(formulas = pulseData$formulas, par = par)
    mean_indexes <-
      sapply(pulseData$conditions, match, names(pulseData$formulas))
    lambdas <- t(t(means[, mean_indexes]) * norm_factors)
    - sum(dnbinom(
      x    = pulseData$count_data,
      mu   = lambdas,
      log  = TRUE,
      size = par$size
    ))
  }
}

ll_norm_factors <- function(pulseData, par) {
  means <- getMeans(formulas = pulseData$formulas, par = par)
  mean_indexes <-
    sapply(pulseData$conditions, match, names(pulseData$formulas))
  lambdas <- means[, mean_indexes]
  norm_indexes <- as.integer(pulseData$fraction)
  function(fraction_factors) {
    fraction_factors <-
      c(1, fraction_factors)[norm_indexes] * pulseData$norm_factors
    norm_lambdas <- t(t(lambdas) * fraction_factors)
    if(any(norm_lambdas<=0)) return(Inf)
    - sum(dnbinom(
      x    = pulseData$count_data,
      mu   = norm_lambdas,
      log  = TRUE,
      size = par$size
    ))
  }
}


ll_dispersion <- function(pulseData, par) {
  norm_factors <- getNormFactors(pulseData, par)
  means <- getMeans(formulas = pulseData$formulas, par = par)
  mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
  lambdas <- t(t(means[, mean_indexes]) * norm_factors)
  function(size) {
    -sum(dnbinom(
        x    = pulseData$count_data,
        mu   = lambdas,
        log  = TRUE,
        size = size
      )
    )
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
  par$fraction_factors <- par$fraction_factors[-1]
  norm_factors <- getNormFactors(pulseData, par)
  means <- getMeans(formulas =  pulseData$formulas, par = par)
  llog <- NULL
  if (!missing(pulseData)) {
    mean_indexes <-
      sapply(pulseData$conditions, match, names(pulseData$formulas))
    lambdas <- t(t(means[, mean_indexes]) * norm_factors)
    llog <- dnbinom(
      x = pulseData$count_data,
      mu = lambdas,
      log = TRUE,
      size = par$size
    )
  }
  list(predictions = lambdas, llog = llog)
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}


evaluateLikelihood <- function(pulseData, par ) {
  objective <- ll_dispersion(pulseData, par)
  objective(par$size)
}

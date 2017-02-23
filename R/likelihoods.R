# extend boundaries to param length
.b <- function(b, par) {
  for (p in names(b)) {
    if (length(b[[p]]) == 1)
      b[[p]] <- rep(b[[p]], length(par[[p]]))
  }
  b
}

#from pryr package
substitute_q <- function(x, env)
{
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

# get matrix for samples
sample_means <- function(evaled_forms, form_indexes, norm_factors){
  mus <- mapply(
    function(i, n){
      m <- do.call(cbind,evaled_forms[i])
      m %*% n
    },
    form_indexes,
    norm_factors, SIMPLIFY=FALSE)
  do.call(cbind, mus)
}

# universal likehood
ll <- function(par, namesToOptimise, pd, singleValue = FALSE) {
  pattern <- par[namesToOptimise]
  if (singleValue)
    pattern <- lapply(pattern, '[[', 1)
  par[namesToOptimise] <- NULL
  function(x, counts) {
    par[namesToOptimise] <- relist(x, pattern)
    evaledForms <- lapply(pd$formulas, eval, envir = par)
    means <- sample_means(evaledForms, pd$formulaIndexes, par$normFactors)
    -sum(dnbinom(counts, mu = means, size = par$size, log = TRUE))
  }
}

# likelihood for norm factors
llnormFactors <- function(par, pd) {
  evaledForms <- lapply(pd$formulas, eval, envir = par)
  function(x, counts) {
    x <- relist(x, par$normFactors)
    means <- sample_means(evaledForms, pd$formulaIndexes, x)
    -sum(dnbinom(counts, mu = means, size = par$size, log = TRUE))
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
  norm_factors <- getNormFactors(pulseData, par)
  means <- getMeans(formulas =  pulseData$formulas, par = par)
  mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
  lambdas <- t(t(means[, mean_indexes]) * norm_factors)
  llog <- dnbinom(
    x = pulseData$count_data,
    mu = lambdas,
    log = TRUE,
    size = par$size
  )
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

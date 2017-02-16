
#' Fit gene-specific parameters
#'
#' @param pulseData PulseData object
#' @param par all needed parameters as a \code{list}
#' @param options \code{list} of options
#'
#' @return fitted parameters as a data.frame ordered as initial parameters
#' @importFrom  stats optim dnbinom
#' @importFrom parallel mclapply
#' 
fitGeneParameters <- function(pulseData, par, options) {
  validateOptions(options)
  param_names <- colnames(par$params)
  lb <- validate(par$params, options$lb$params)
  ub <- validate(par$params, options$ub$params)
  objective <- ll_gene(pulseData, par)
  new_params <- list()
  new_params <- mclapply(
    X = seq_len(dim(par$params)[1]),
    FUN = function(i) {
      olds <- par$params[i,, drop = FALSE]
      optim(
        olds,
        objective,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(parscale = olds),
        counts = pulseData$count_data[i,],
        known = par$known[i,, drop = FALSE]
      )$par
    }
    ,mc.cores = options$cores
  )
  new_params <- do.call(rbind, new_params)
  rownames(new_params) <- rownames(par$params)
  as.data.frame(new_params)
}

fitSharedParameters <- function(pulseData, par, options) {
  shared_objective <- ll_shared_params(pulseData, par)
  shared_names <- names(par$shared)
  lb <- validate(par$shared, options$lb$shared)
  ub <- validate(par$shared, options$ub$shared)
  shared_params <- optim(
    unlist(par$shared),
    shared_objective,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub
  )$par
  names(shared_params) <- shared_names
  as.list(shared_params)
}

#' Fit the dispersion parameter
#'
#' @param pulseData PulseData object
#' @param par all needed parameters as a \code{list}
#' @param options \code{list} of options
#'
#' @return the MLE of the \code{size} for \code{\link{dnbinom}}
#' @importFrom  stats optimise
#'
fitDispersion <- function(pulseData, par, options) {
  dispersion_objective <- ll_dispersion(pulseData, par)
  interval <- c(options$lb$size, options$ub$size)
  size <- optimise(dispersion_objective,
                   interval = interval)$minimum
  size
}

#' Fit fraction normalisation coefficients
#'
#' @inheritParams fitDispersion
#' @return vector of normalisation factors; \code{c(1,norm_factors)}
#'   corresponds to the fractions in \code{pulseData$fractions}
#'
#' @importFrom  stats optimise
fitFractions <- function(pulseData, par, options){
  objective <- ll_norm_factors(pulseData, par)
  fraction_names <- names(par$fraction_factors)
  lb <- validate(par$fraction_factors, options$lb$fraction_factors)
  ub <- validate(par$fraction_factors, options$ub$fraction_factors)
  fraction_factors <- optim(
    unlist(par$fraction_factors)[-1],
    objective,
    method = "L-BFGS-B",
    lower  = options$lb$fraction_factors,
    upper  = options$ub$fraction_factors
  )$par
  c(1,fraction_factors)
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

#' Fit the model by MLE
#'
#' @param pulseData PulseData object
#' @param par initial guess for parameters as a \code{list}.  
#' @param options \code{list} of options
#'
#' @return a list with the fitted parameters in the same form as
#' the initial guess `par`
#'     
#' @details 
#'    In the initial guess \code{par},
#'    the following list items need to be provided:
#'    \itemize{
#'    \item{params: }{list of gene-specific parameters}
#'    \item{shared: }{list of shared parameters}
#'    \item{fraction_factors (if relevant): }{a vector of fraction factors
#'    in the order of \code{levels(pulseDatad$fractions)}} 
#'    \item{size (optional): }{an initial guess for the size parameter of the
#'    negative binomial distribution, see \code{dnbinom}}
#'    }
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options) {
  param_names <- names(par$params)
  log2screen(options, cat("\n"))
  rel_err <- Inf
  shared_rel_err <- ifelse(is.null(par$shared), 0, Inf)
  fraction_rel_err <- ifelse(is.null(par$fraction_factors), 0, Inf)
  while (rel_err > options$tolerance$params ||
         shared_rel_err > options$tolerance$shared ||
         fraction_rel_err > options$tolerance$fraction) {
    # Fit shared params
    if (!is.null(par$shared_params)) {
      shared_params <- fitSharedParameters(pulseData, par, options)
      shared_rel_err <- getMaxRelDifference(shared_params, par$shared)
      par$shared <- shared_params
    }
    # Fit params for every genes individually
    params <- fitGeneParameters(pulseData, par, options)
    rel_err <- getMaxRelDifference(params, par$params)
    par$params <- params
    if (!is.null(par$fraction_factors)) {
      res <- fitFractions(pulseData, par, options)
      fraction_rel_err <- getMaxRelDifference(res, par$fraction_factors)
      par$fraction_factors <- res
    }
    par$size <- fitDispersion(pulseData, par, options)
    str <- format(c(rel_err, shared_rel_err, fraction_rel_err),
                  digits = 2,
                  width = 6)
    log2screen(options, cat(
      paste0(
        "Max Rel.err. in [params: ", str[1],
        "]  [shared: ", str[2], 
        "]  [fractions: ", str[3],
        "]    \r"
      )
    ))
  }
  ## fit gene specific final parameters
  par$individual_params <- fitGeneParameters(pulseData, par, options)
  names(par$fraction_factors)  <- levels(pulseData$fraction)
  par
}

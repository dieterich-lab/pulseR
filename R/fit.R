
.defaultParams <- function() {
  list(
    rel_tol = 1e-3,
    shared_rel_tol = 1e-2,
    fraction_rel_err = 1e-3,
    verbose = "silent",
    update_inital_parameters = FALSE,
    cores = 1,
    lower_boundary_size = 10,
    upper_boundary_size = 1e10
  )
}
 
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
fitIndividualParameters <- function(pulseData, par, options) {
  param_names <- colnames(par$individual_params)
  objective <- ll_gene(pulseData, par)
  new_params <- list()
  new_params <- mclapply(
    X = seq_len(dim(par$individual_params)[1]),
    FUN = function(i) {
      olds <- par$individual_params[i,]
      optim(
        olds,
        objective,
        method = "L-BFGS-B",
        lower = options$lower_boundary,
        upper = options$upper_boundary,
        control = list(parscale = olds),
        counts = pulseData$count_data[i,],
        known=par$known[i,]
      )$par
    }
    ,mc.cores = options$cores
  )
  new_params <- do.call(rbind, new_params)
  rownames(new_params) <- rownames(par$individual_params)
  as.data.frame(new_params)
}

fitSharedParameters <- function(pulseData, par, options) {
  shared_objective <- ll_shared_params(pulseData, par)
  shared_params <- optim(
    unlist(par$shared_params),
    shared_objective,
    method = "L-BFGS-B",
    lower = options$lower_boundary_shared,
    upper = options$upper_boundary_shared
  )$par
  names(shared_params) <- names(par$shared_params)
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
  size <- optimise(dispersion_objective,
                   interval = unlist(options[c("lower_boundary_size",
                                               "upper_boundary_size")]))$minimum
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
  fraction_factors <- optim(
    unlist(par$fraction_factors),
    objective,
    method = "L-BFGS-B",
    lower = options$lower_boundary_fraction,
    upper = options$upper_boundary_fraction
  )$par
  fraction_factors
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

#' Fit the model by MLE
#'
#' @param pulseData PulseData object
#' @param par all needed parameters as a \code{list}
#' @param options \code{list} of options
#'
#' @return a list `l` with fitted parameters `l$par` and
#'     formulas `l$formulas` 
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options = list()) {
  opts <- .defaultParams()
  param_names <- names(par$individual_params)
  opts$parscales <- mapply(max,
    abs(options$upper_boundary),
    abs(options$lower_boundary))
  rel_err <- 10 * opts$rel_tol
  if (is.null(par$shared_params)) {
    shared_rel_err <- 0
  } else {
    shared_params <- as.list(par$shared_params)
    shared_rel_err <- 10 * opts$shared_rel_tol
  }
  if (is.null(par$fraction_factors)){
    fraction_rel_err <- 0 
  } else {
    fraction_rel_err <- 10*opts$fraction_rel_err
  }
  opts[names(options)] <- options
  while (rel_err > opts$rel_tol ||
         shared_rel_err > opts$shared_rel_tol ||
         fraction_rel_err > opts$fraction_rel_err) {
    # Fit shared params
    if (!is.null(par$shared_params)) {
      log2screen(opts, "Fitting shared params\n")
      shared_params <- fitSharedParameters(pulseData, par, opts)
      shared_rel_err <- getMaxRelDifference(shared_params, par$shared_params)
      par$shared_params <- shared_params
      log2screen(opts, "Shared params\n")
      log2screen(opts, toString(par$shared_params), "\n")
    }
    # Fit params for every genes individually
    log2screen(opts, "Fitting gene-specific params\n")
    params <- fitIndividualParameters(pulseData, par, opts)
    rel_err <- getMaxRelDifference(params, par$individual_params)
    par$individual_params <- params
    if (!is.null(par$fraction_factors)) {
      log2screen(opts, "Fitting fraction coefficients\n")
      res <- fitFractions(pulseData, par, opts)
      fraction_rel_err <- getMaxRelDifference(res, par$fraction_factors)
      par$fraction_factors <- res
      log2screen(opts, "fraction factors \n")
      log2screen(opts, toString(par$fraction_factors), "\n")
    }
    par$size <- fitDispersion(pulseData, par, opts)
  }
  ## fit gene specific final parameters
  par$individual_params <- fitIndividualParameters(pulseData, par, opts)
  par$fraction_factors <- c(1, par$fraction_factors)
  names(par$fraction_factors)  <- levels(pulseData$fraction)
  list(par = par, formulas = pulseData$formulas)
}

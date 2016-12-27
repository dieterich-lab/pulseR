
defaultParams <- function() {
  list(
    rel_tol = 1e-2,
    shared_rel_tol = 1e-2,
    verbose = "silent",
    update_inital_parameters = FALSE,
    cores = 1,
    lower_boundary_size = 10,
    upper_boundary_size = 1e10
  )
}
 
## params are a data.frame with ids in rownames
# order is the same as in parameters
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
        control = list(parscale = options$parscales),
        counts = pulseData$count_data[i,]
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

fitDispersion <- function(pulseData, par, options) {
  dispersion_objective <- ll_dispersion(pulseData, par)
  size <- optimise(dispersion_objective,
                   interval = unlist(options[c("lower_boundary_size",
                                               "upper_boundary_size")]))$minimum
  size
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

# options is a list with records
# - individual_rel_err
# - shared_rel_tol
fitModel <- function(pulseData, par, options = list()) {
  opts <- defaultParams()
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
  opts[names(options)] <- options
  while (rel_err > opts$rel_tol ||
         shared_rel_err > opts$shared_rel_tol) {
    # Fit params for every genes individually
    log2screen(opts, "Fitting gene-specific params\n")
    params <- fitIndividualParameters(pulseData, par, opts)
    rel_err <- getMaxRelDifference(params, par$individual_params)
    par$individual_params <- params
    # Fit shared params
    if (!is.null(par$shared_params)) {
      log2screen(opts, "Fitting shared params\n")
      shared_params <- fitSharedParameters(pulseData, par, opts)
      shared_rel_err <- getMaxRelDifference(shared_params, par$shared_params)
      par$shared_params <- shared_params
      log2screen(opts, "Shared params\n")
      log2screen(opts, toString(par$shared_params), "\n")
    }
    par$size <- fitDispersion(pulseData, par, opts)
  }
  list(par=par, formulas = pulseData$formulas)
}


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
fitIndividualParameters <- function(pulseData,
                                    old_params,
                                    shared_params,
                                    options,
                                    size) {
  param_names <- colnames(old_params)
  objective <- ll_gene(
    pulseData = pulseData,
    param_names = param_names,
    size = size,
    shared_params = shared_params
  )
  new_params <- list()
  new_params <- mclapply(
    X = seq_len(dim(old_params)[1]),
    FUN = function(i) {
      olds <- old_params[i,]
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
  rownames(new_params) <- rownames(old_params)
  as.data.frame(new_params)
}

fitSharedParameters <- function(old_shared_params,
                                pulseData,
                                individual_params,
                                options,
                                size) {
  shared_objective <- ll_shared_params(
    pulseData,
    individual_params = individual_params,
    shared_param_names =  names(old_shared_params),
    size =  size
  )
  shared_params <- optim(
    unlist(old_shared_params),
    shared_objective,
    method = "L-BFGS-B",
    lower = options$lower_boundary_shared,
    upper = options$upper_boundary_shared
  )$par
  names(shared_params) <- names(old_shared_params)
  as.list(shared_params)
}

fitDispersion <- function(shared_params,
                          pulseData,
                          individual_params,
                          options,
                          size) {
  dispersion_objective <- ll_dispersion(
    pulseData = pulseData,
    individual_params =  individual_params,
    shared_params =  shared_params
  )
  size <- optimise(dispersion_objective,
                   interval = unlist(options[c("lower_boundary_size",
                                               "upper_boundary_size")]))$minimum
  size
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

# options is a list with records
# - individual_rel_err
# - shared_rel_tol
fitModel <- function(pulseData,
                     params,
                     shared_params = NULL,
                     size=100,
                     options = list()) {
  require(parallel)
  opts <- defaultParams()
  param_names <- names(params)
  opts$parscales <- mapply(max,
                           abs(options$upper_boundary),
                           abs(options$lower_boundary))
  rel_err <- 10 * opts$rel_tol
  if (is.null(shared_params)) {
    shared_rel_err <- 0
  } else {
    shared_params <- as.list(shared_params)
    shared_rel_err <- 10 * opts$shared_rel_tol
  }
  opts[names(options)] <- options
  while (rel_err > opts$rel_tol ||
         shared_rel_err > opts$shared_rel_tol) {
    # Fit params for every genes individually
    log2screen(opts, "Fitting gene-specific params\n")
    old_params <- params
    params <- fitIndividualParameters(
      old_params = old_params,
      pulseData = pulseData,
      shared_params = shared_params,
      options = opts,
      size = size
    )
    rel_err <- getMaxRelDifference(params, old_params)
    # Fit shared params
    if (!is.null(shared_params)) {
      log2screen(opts, "Fitting shared params\n")
      old_shared_params <- shared_params
      shared_params <- fitSharedParameters(
        old_shared_params = shared_params,
        pulseData = pulseData,
        individual_params = params,
        options           = opts,
        size              = size
      )
      shared_rel_err <-
        getMaxRelDifference(shared_params, old_shared_params)
      log2screen(opts, "Shared params\n")
      log2screen(opts, toString(shared_params), "\n")
    }
    size <- fitDispersion(
      shared_params = shared_params,
      pulseData = pulseData,
      individual_params = params,
      options = opts,
      size = size
    )
  }
  list(
    individual_params = params,
    shared_params = shared_params,
    size = size, 
    formulas = pulseData$formulas
  )
}
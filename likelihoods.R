library(parallel)

#from pryr package
substitute_q <- function (x, env)
{
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

MeanFormulas <- function(...) {
  eval(substitute(alist(...)))
}

makeVector <- function(forms) {
  string <- paste("c(",
                  paste(mapply(
                    paste, names(forms), as.character(forms), sep = "="
                  ),
                  collapse = ","),
                  ")")
  parse(text = string)[[1]]
}

contaminate <- function(formulas,
                        target_condition,
                        contaminant_condition,
                        coef_name) {
  target_condition <- deparse(substitute(target_condition))
  contaminant_condition <-
    deparse(substitute(contaminant_condition))
  f1 <- deparse(formulas[[target_condition]])
  f2 <- deparse(formulas[[contaminant_condition]])
  e <- paste("(1-", coef_name, ")*(", f1, ")+", coef_name, "*(", f2, ")")
  parse(text = e)[[1]]
}

constructFormulas <- function(formulas, conditions) {
  result <- lapply(rownames(conditions), function(x) {
    sampleCondition <- conditions[x, 1L]
    substitute_q(formulas[[sampleCondition]],
                 conditions[x, -1, drop = FALSE])
  })
  names(result) <- rownames(conditions)
  result
}

ll_gene <- function(conditions,
                    formulas,
                    param_names,
                    size,
                    norm_factors = rep(1, dim(counts)[2]),
                    shared_params = NULL) {
  mean_indexes <- sapply(conditions, match, names(formulas))
  if (!is.null(shared_params) && !anyNA(formulas))
    formulas <- lapply(formulas, substitute_q, shared_params)
  means_vector <- makeVector(formulas)
  funquote <- function(params, counts) {
    names(params) <- param_names
    mus <- eval(means_vector, as.list(params))
    lambdas <- norm_factors * mus[mean_indexes] + 1e-10
    -sum(dnbinom(
      x    = counts,
      mu   = lambdas, 
      log  = TRUE,
      size = size
    ))
  }
  funquote
}

ll_shared_params <- function(count_data,
                             conditions,
                             norm_factors,
                             formulas,
                             individual_params,
                             shared_param_names,
                             size) {
  function(shared_params) {
  names(shared_params) <- shared_param_names
    means <- getMeans (shared_params,
                       formulas,
                       individual_params)
    mean_indexes <- sapply(conditions, match, names(formulas))
    lambdas <- means[, mean_indexes]
    - sum(dnbinom(
      x    = count_data,
      mu   = lambdas * norm_factors,
      log  = TRUE,
      size = size
    ))
  }
}

getMeans <- function(shared_params,
                     formulas,
                     individual_params) {
  shared_params <- as.list(shared_params)
  means <- lapply(formulas, function(x) {
    eval(substitute_q(x, shared_params),
         envir = as.list(individual_params))
  })
  means <- do.call(cbind, means) + 1e-10
  means
}

ll_dispersion <- function(count_data,
                          conditions,
                          norm_factors,
                          formulas,
                          individual_params,
                          shared_params) {
  function(size) {
    means <- getMeans(shared_params,
                      formulas,
                      individual_params)
    mean_indexes <- sapply(conditions, match, names(formulas))
    lambdas <- means[, mean_indexes]
    - sum(dnbinom(
      x = count_data,
      mu = lambdas * norm_factors,
      log = TRUE,
      size = size
    ))
  }
}

predict.expression <- function(count_data,
                               conditions,
                               formulas,
                               fit) {
  means <- getMeans(fit$shared_params,
                    formulas,
                    fit$individual_params)
  mean_indexes <- sapply(conditions, match, names(formulas))
  lambdas <- means[, mean_indexes]
  llog <- -(dnbinom(
    x = count_data,
    mu = lambdas * fit$norm_factors,
    log = TRUE,
    size = fit$size
  ))
  list(preditions = lambdas, llog = llog)
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}

## params are a data.frame with ids in rownames
# order is the same as in parameters
fitIndividualParameters <- function(old_params,
                                    count_data,
                                    conditions,
                                    formulas,
                                    shared_params,
                                    norm_factors,
                                    options,
                                    size) {
  param_names <- colnames(old_params)
  objective <- ll_gene(
    conditions = conditions$condition,
    formulas = formulas,
    param_names = param_names,
    size = size,
    norm_factors = norm_factors,
    shared_params = shared_params
  )
  new_params <- list()
  new_params <- mclapply(
    X = seq_len(dim(old_params)[1]),
    FUN = function(i) {
      olds <- old_params[i, ]
      optim(
        olds,
        objective,
        method = "L-BFGS-B",
        lower = options$lower_boundary,
        upper = options$upper_boundary,
        control = list(parscale = olds),
        counts = count_data[i, ]
      )$par
    },
    mc.cores = options$cores
  )
  new_params <- do.call(rbind, new_params)
  rownames(new_params) <- rownames(old_params)
  as.data.frame(new_params)
}

fitSharedParameters <- function(old_shared_params,
                                count_data,
                                conditions,
                                formulas,
                                individual_params,
                                norm_factors,
                                options,
                                size) {
  shared_objective <- ll_shared_params(
    count_data = count_data,
    conditions = conditions$condition,
    norm_factors = norm_factors,
    formulas = formulas,
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

evaluateLikelihood <- function(shared_params,
                               count_data,
                               conditions,
                               formulas,
                               individual_params,
                               norm_factors,
                               size) {
  shared_objective <- ll_shared_params(
    count_data = count_data,
    conditions = conditions$condition,
    norm_factors = norm_factors,
    formulas = formulas,
    individual_params = individual_params,
    shared_param_names =  names(shared_params),
    size =  size
  )
  shared_objective(unlist(shared_params))
}

fitDispersion <- function(shared_params,
                          count_data,
                          conditions,
                          norm_factors,
                          formulas,
                          individual_params,
                          options,
                          size) {
  dispersion_objective <- ll_dispersion(
    count_data = count_data,
    conditions = conditions,
    norm_factors = norm_factors,
    formulas =  formulas,
    individual_params =  individual_params,
    shared_params =  shared_params
  )
  size <- optimise(dispersion_objective,
                   interval = unlist(options[c("lower_boundary_size",
                                               "upper_boundary_size")]))$minimum
  size
}

# options is a list with records
# - individual_rel_err
# - shared_rel_tol
fitModel <- function(count_data,
                     conditions,
                     norm_factors,
                     formulas,
                     params,
                     shared_params = NULL,
                     options = list()) {
  require(parallel)
  param_names <- names(params)
  opts <- list(
    rel_tol = 1e-2,
    shared_rel_tol = 1e-2,
    verbose = "silent",
    update_inital_parameters = FALSE,
    cores = 1,
    lower_boundary_size = 10,
    upper_boundary_size = 1e10
  )
  rel_err <- 10 * opts$rel_tol
  if (is.null(shared_params)) {
    shared_rel_err <- 0
  } else {
    shared_rel_err <- 10 * opts$shared_rel_tol
  }
  size <- 1e2
  opts[names(options)] <- options
  while (rel_err > opts$rel_tol ||
         shared_rel_err > opts$shared_rel_tol) {
    # Fit params for every genes individually
    log2screen(opts, "Fitting gene-specific params")
    old_params <- params
    params <- fitIndividualParameters(
      old_params = old_params,
      count_data = count_data,
      conditions = conditions,
      formulas = formulas,
      shared_params = shared_params,
      norm_factors = norm_factors,
      options = opts,
      size = size
    )
    rel_err <- getMaxRelDifference(params, old_params)
    # Fit shared params
    if (!is.null(shared_params)) {
      log2screen(opts, rep(" ", 100, "\r"))
      log2screen(opts, "Fitting shared params")
      old_shared_params <- shared_params
      shared_params <- fitSharedParameters(
        old_shared_params = shared_params,
        count_data        = count_data,
        conditions        = conditions,
        formulas          = formulas,
        individual_params = params,
        norm_factors      = norm_factors,
        options           = opts,
        size              = size
      )
      shared_rel_err <- getMaxRelDifference(shared_params, old_shared_params)
      log2screen(opts, "Shared params\n")
      log2screen(opts, toString(shared_params), "\n")
    }
    size <- fitDispersion(
      shared_params = shared_params,
      count_data = count_data,
      conditions = conditions,
      formulas = formulas,
      individual_params = params,
      norm_factors = norm_factors,
      options = opts,
      size = size
    )
  }
  list(
    individual_params = params,
    shared_params = shared_params,
    size = size,
    norm_factors = norm_factors
  )
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

findDeseqFactorsSingle <- function(count_data, spikeins)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    exp(median(log(x) - loggeomeans, na.rm = TRUE))
  })
  deseqFactors
}

findDeseqFactors <- function(count_data, conditions, spikeins) {
  lapply(split(rownames(conditions), conditions$condition),
         function(samples) {
           findDeseqFactorsSingle(count_data[, samples], spikeins)
         })
}
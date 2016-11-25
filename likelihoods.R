library(parallel)

#from pryr package
substitute_q <- function (x, env)
{
  stopifnot(is.language(x))
  env <- to_env(env)
  call <- substitute(substitute(x, env), list(x = x))
  eval(call)
}

# from pryr package
to_env <-  function (x, quiet = FALSE)
{
  if (is.environment(x)) {
    x
  }
  else if (is.list(x)) {
    list2env(x)
  }
  else if (is.function(x)) {
    environment(x)
  }
  else if (length(x) == 1 && is.character(x)) {
    if (!quiet)
      message("Using environment ", x)
    as.environment(x)
  }
  else if (length(x) == 1 && is.numeric(x) && x > 0) {
    if (!quiet)
      message("Using environment ", search()[x])
    as.environment(x)
  }
  else {
    stop("Input can not be coerced to an environment", call. = FALSE)
  }
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
  parse(text = string)
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
  e <-
    paste("(1-", coef_name, ")*(", f1, ")+", coef_name, "*(", f2, ")")
  parse(text = e)[[1]]
}


ll_gene <- function(count_data,
                    formulas,
                    param_names,
                    size,
                    shared_params = NULL) {
  conditions <- count_data$condition
  counts <- count_data$count
  norm_factors <- count_data$norm_factor
  stopifnot(!is.null(conditions),
            !is.null(norm_factors),
            !is.null(counts))
  mean_indexes <- sapply(conditions, match, names(formulas))
  if (!is.null(shared_params) && !anyNA(formulas))
    formulas <- lapply(formulas, substitute_q, shared_params)
  means_vector <- makeVector(formulas)
  funquote <- substitute(function(params) {
    names(params) <- param_names
    mus <- eval(means_vector, as.list(params))
    lambdas <- mus[mean_indexes] + 1e-10
    - sum(dnbinom(
      counts,
      mu = lambdas * norm_factors,
      log = TRUE,
      size = size
    ))
  }, parent.frame())
  eval(funquote)
}

getMeansEstimatingFunction <- function(count_data,
                                       individual_params,
                                       forms,
                                       shared_param_names) {
  gene_number <- length(individual_params$id)
  individual_params$id <- as.factor(individual_params$id)
  forms_vector <- makeVector(forms)[[1]]
  condition_indexes <-
    sapply(count_data$condition, match, names(forms))
  means_indexes <-
    cbind(as.numeric(count_data$id), condition_indexes)
  means_for_genes <- list()
  for (i in 1:gene_number) {
    means_for_genes[[i]] <-
      substitute_q(forms_vector, env = individual_params[i, ])
  }
  means_matrix <-
    matrix(0, ncol = length(forms), nrow = gene_number)
  function(shared_params) {
    names(shared_params) <- shared_param_names
    means_matrix <- do.call(rbind, lapply(means_for_genes, eval,
                                          envir = as.list(shared_params)))
    means_matrix[means_indexes]
  }
}

ll_shared_params <- function(count_data,
                             forms,
                             individual_params,
                             shared_param_names,
                             size) {
  estimateMeans <-
    getMeansEstimatingFunction(count_data, individual_params,
                               forms, shared_param_names)
  norm_factors <- count_data$norm_factor
  function(shared_params) {
    lambdas <- estimateMeans(shared_params) + 1e-10
    -sum(dnbinom(
      x = count_data$count,
      mu = lambdas * norm_factors,
      log = TRUE,
      size = size
    ))
  }
}

ll_dispersion <- function(count_data,
                          forms,
                          individual_params,
                          shared_params,
                          size) {
  estimateMeans <-
    getMeansEstimatingFunction(count_data, individual_params,
                               forms, names(shared_params))
  lambdas <- estimateMeans(shared_params) + 1e-10
  lambdas <- lambdas * count_data$norm_factor
  function(size) {
    -sum(dnbinom(
      x = count_data$count,
      mu = lambdas,
      log = TRUE,
      size = size
    ))
  }
}


predict.expression <- function(count_data, model, forms) {
  estimateMeans <- getMeansEstimatingFunction(count_data,
                                              model$individual_params,
                                              forms,
                                              names(model$shared_params))
  lambdas <- estimateMeans(model$shared_params)
  list(prediction = lambdas,
       logL = dpois(count_data$count, (lambdas + 1e-10), log = TRUE))
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}

fitIndividualParameters <- function(old_params,
                                    splitted_data,
                                    formulas,
                                    shared_params,
                                    options,
                                    size) {
  param_names <- setdiff(names(old_params), "id")
  objectives <- lapply(
    splitted_data,
    ll_gene,
    formulas = formulas,
    param_names = param_names,
    size = size,
    shared_params = shared_params
  )
  ids <- old_params$id
  old_params <- split(old_params[, param_names],
                      old_params$id)[names(objectives)]
  new_params <- list()
  new_params <- mcmapply(
    function(obj, olds) {
      optim(
        unlist(olds),
        obj,
        method = "L-BFGS-B",
        lower = options$lower_boundary,
        upper = options$upper_boundary
      )$par
    },
    objectives,
    old_params,
    SIMPLIFY = FALSE,
    mc.cores = options$cores
  )
  new_params <- as.data.frame(do.call(rbind, new_params))
  names(new_params) <- param_names
  new_params$id <- ids
  new_params
}

fitSharedParameters <- function(old_shared_params,
                                count_data,
                                formulas,
                                individual_params,
                                options,
                                size) {
  shared_objective <- ll_shared_params(count_data,
                                       formulas,
                                       individual_params,
                                       names(old_shared_params),
                                       size)
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
                               individual_params,
                               count_data,
                               formulas) {
  conditions <- names(formulas)
  count_data <- count_data[count_data$condition %in% conditions, ]
  shared_objective <- ll_shared_params(count_data,
                                       formulas,
                                       individual_params,
                                       names(shared_params))
  shared_objective(unlist(shared_params))
}

fitDispersion <- function(shared_params,
                          count_data,
                          formulas,
                          individual_params,
                          options,
                          size) {
  dispersion_objective <- ll_dispersion(count_data,
                                        formulas,
                                        individual_params,
                                        shared_params,
                                        size)
  size <- optimise(dispersion_objective,
                   interval = unlist(options[c("lower_boundary_size",
                                               "upper_boundary_size")]))$minimum
  size
}

# options is a list with records
# - individual_rel_err
# - shared_rel_tol
fitModel <- function(count_data,
                     formulas,
                     individual_params,
                     shared_params = NULL,
                     options = list()) {
  require(parallel)
  conditions <- names(formulas)
  count_data <- count_data[count_data$condition %in% conditions, ]
  splitted_data <- split(count_data, count_data$id)
  param_names <- setdiff(names(individual_params), "id")
  shared_param_names <- names(shared_params)
  opts <- list(
    individual_rel_tol = 1e-2,
    shared_rel_tol = 1e-2,
    verbose = "silent",
    update_inital_parameters = FALSE,
    cores = 1,
    lower_boundary_size = 10,
    upper_boundary_size = 1e10
  )
  individual_rel_err <- 10 * opts$individual_rel_tol
  if (is.null(shared_params)) {
    shared_rel_err <- 0
  } else {
    shared_rel_err <- 10 * opts$shared_rel_tol
  }
  size <- 1000
  opts[names(options)] <- options
  while (individual_rel_err > opts$individual_rel_tol ||
         shared_rel_err > opts$shared_rel_tol) {
    shared_params <- as.list(shared_params)
    # Fit params for every genes individually
    log2screen(opts, "Fitting gene-specific params")
    old_params <- individual_params
    individual_params <-
      fitIndividualParameters(
        old_params = old_params,
        splitted_data = splitted_data,
        formulas = formulas,
        shared_params = shared_params,
        options = opts,
        size = size
      )
    individual_rel_err <-
      max(abs(1 - individual_params[, param_names] / old_params[, param_names]))
    # Fit shared params
    if (!is.null(shared_params)) {
      log2screen(opts, rep(" ", 100, "\r"))
      log2screen(opts, "Fitting shared params")
      old_shared_params <- shared_params
      shared_params <-
        fitSharedParameters(
          old_shared_params = shared_params,
          count_data        = count_data,
          formulas          = formulas,
          individual_params = individual_params,
          options           = opts,
          size              = size
        )
      shared_rel_err <- max(abs(1 - unlist(shared_params) / unlist(old_shared_params)))
      log2screen(opts, "Shared params\n")
      log2screen(opts, toString(shared_params), "\n")
    }
    size <- fitDispersion(
      shared_params = shared_params,
      count_data = count_data,
      formulas = formulas,
      individual_params = individual_params,
      options = opts,
      size = size
    )
  }
  list(
    individual_params = individual_params,
    shared_params = shared_params,
    size = size
  )
}

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

ll_gene <- function(pulseData, par) {
  mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
  formulas <- pulseData$formulas
  norm_factors <- pulseData$norm_factors
  if(!is.null(par$fraction_factors)){
    norm_factors <- 
      norm_factors * par$fraction_factors[as.integer(pulseData$conditions$fraction)]
  }
  if (!is.null(par$shared_params))
    formulas <- lapply(formulas, substitute_q, par$shared_params)
  means_vector <-  makeVector(formulas)
  param_names <- names(par$individual_params)
  funquote <- function(params, counts) {
    names(params) <- param_names
    mus <- eval(means_vector, as.list(params))
    lambdas <-  mus[mean_indexes] 
    -sum(dnbinom(
      x    = counts,
      mu   = lambdas * norm_factors,
      log  = TRUE,
      size = par$size
    ))
  }
  funquote
}

ll_shared_params <- function(pulseData, par) {
  norm_factors <- pulseData$norm_factors 
  if(!is.null(par$fraction_factors)){
    norm_factors <- 
      norm_factors * par$fraction_factors[as.integer(pulseData$conditions$fraction)]
  }
  function(shared_params) {
    names(shared_params) <- shared_param_names
    means <- getMeans(shared_params,
                      pulseData$formulas,
                      par$individual_params)
    mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
    lambdas <- means[, mean_indexes]
    - sum(
      dnbinom(
        x    = pulseData$count_data,
        mu   = lambdas * norm_factors,
        log  = TRUE,
        size = par$size
      )
    )
  }
}

ll_norm_factors <- function(pulseData, par) {
  means <- getMeans(par$shared_params,
    pulseData$formulas,
    par$individual_params)
  mean_indexes <- sapply(pulseData$conditions, match, names(pulseData$formulas))
  norm_lambdas <- means[, mean_indexes] * pulseData$norm_factors
  norm_indexes <- as.integer(conditions$fraction)
  function(fraction_factors) {
    fraction_factors <- c(1,fraction_factors)
    - sum(
      dnbinom(
        x    = pulseData$count_data,
        mu   = norm_lambdas * fraction_factors[norm_indexes],
        log  = TRUE,
        size = par$size
      )
    )
  }
}

getMeans <- function(shared_params, formulas, individual_params) {
  shared_params <- as.list(shared_params)
  means <- lapply(formulas, function(x) {
    eval(substitute_q(x, shared_params),
         envir = as.list(individual_params))
  })
  means <- do.call(cbind, means) 
  means
}

ll_dispersion <- function(pulseData,
                          individual_params,
                          shared_params) {
  function(size) {
    means <- getMeans(shared_params,
                      pulseData$formulas,
                      individual_params)
    mean_indexes <-
      sapply(pulseData$conditions, match, names(pulseData$formulas))
    lambdas <- means[, mean_indexes]
    -sum(dnbinom(
        x    = pulseData$count_data,
        mu   = lambdas * pulseData$norm_factors,
        log  = TRUE,
        size = size
      )
    )
  }
}

predict.expression <- function(fit, pulseData) {
  means <- getMeans(fit$shared_params,
                    fit$formulas,
                    fit$individual_params)
  llog <- NULL
  if (!missing(pulseData)) {
    mean_indexes <-
      sapply(pulseData$conditions[, 1], match, names(pulseData$formulas))
    means <- means[, mean_indexes] * pulseData$norm_factors
    colnames(means) <- colnames(pulseData$count_data)
    llog <- dnbinom(
      x = pulseData$count_data,
      mu = means,
      log = TRUE,
      size = fit$size
    )
  }
  list(predictions = means, llog = llog)
}

log2screen <- function(options, ...) {
  if (options$verbose == "verbose")
    cat(...)
}


evaluateLikelihood <- function(shared_params,
                               pulseData,
                               individual_params,
                               size) {
  shared_objective <- ll_shared_params(
    pulseData,
    individual_params = individual_params,
    shared_param_names =  names(shared_params),
    size =  size
  )
  shared_objective(unlist(shared_params))
}

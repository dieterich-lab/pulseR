source("likelihoods.R")

generateTestDataSingle <- function(forms,
                                   individual_params,
                                   shared_params,
                                   conditions = rep(names(forms), n),
                                   n = 1,
                                   size = 100) {
  means <- sapply(forms, eval,
                  c(as.list(individual_params), as.list(shared_params)))
  indexes <- match(conditions, names(forms))
  counts <- rnbinom(n    = length(conditions),
                    mu   = means[indexes],
                    size = size)
  counts
}

generateTestData <- function(n, replicates) {
  conditions <- list()
  conditions$sample <- replicate(length(forms) * replicates,
                                 paste0(letters[sample(25, 10)], collapse = ""))
  conditions$condition <- names(forms)
  conditions <- as.data.frame(conditions, stringAsFactors = FALSE)
  genes <- replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- paste0("ENS00000", genes)
  p <- data.frame(
    id = genes,
    mu_n = runif(n, 1e2, 50000),
    mu_h = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8),
    a_h  = runif(n, .05, .8)
  )
  alphas <- list(
    alpha_chase = 1,
    alpha_lab = 2,
    beta_chase = 1,
    beta_lab = 2
  )
  size <- 1e2
  d <- lapply(seq_along(p$id), function(i) {
    generateTestDataSingle(forms, p[i, -1], alphas, n = replicates)
  })
  data <- do.call(rbind, d)
  rownames(data) <- genes
  colnames(data) <- conditions$sample
  rownames(conditions) <- conditions$sample
  list(
    data = data,
    conditions = conditions[, -1, drop = FALSE],
    params = p,
    size = size,
    shared_params = alphas
  )
}

forms <- MeanFormulas(
  total_Norm         = mu_n,
  total_Hypox        = mu_n * a_h + mu_h,
  flow_lab_Norm      = alpha_lab * (mu_n * a_n),
  biotin_lab_Norm    = beta_lab * mu_n * (1 - a_n),
  flow_lab_Hypox     = alpha_lab * (mu_n * a_h),
  biotin_lab_Hypox   = beta_lab * mu_h,
  flow_chase_Norm    = alpha_chase * mu_n * (1 - (1 - a_n) * a_n),
  biotin_chase_Norm  = beta_chase * mu_n * (1 - a_n) * a_n,
  flow_chase_Hypox   = alpha_chase * (mu_h + mu_n * a_n * a_h),
  biotin_chase_Hypox = beta_chase * (mu_n * (1 - a_n) * a_h)
)

guess_params <- function(data, params) {
  mean_expression <-
    unlist(lapply(data, function(x)
      mean(x$count[x$condition == "total_Norm"])))
  guess <- data.frame(
    id = params$id,
    mu_n = mean_expression,
    mu_h = mean_expression,
    a_n = .1,
    a_h = .1
  )
  guess
}

testIndividualGeneParams <- function(n = 2, replicates = 2) {
  g <- generateTestData(n = n, r = replicates)
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    cores = 2
  )
  g$data <- g$data[order(g$data$id),]
  g$params <- g$params[order(g$params$id),]
  data <- split(g$data, g$data$id)
  guess <- guess_params(data, g$params)
  estimation <- fitIndividualParameters(guess, data, forms, g$shared_params,
                                        options, g$size)
  errors <-
    abs(1 - estimation[, which(names(estimation) != "id"), drop = FALSE] /
          g$params[, which(names(g$params) != "id"), drop = FALSE])
  unlist(errors)
}

testSharedParams <- function(n = 2, replicates = 2) {
  g <- generateTestData(n = n, r = replicates)
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    cores = 2
  )
  g$data <- g$data[order(g$data$id),]
  g$params <- g$params[order(g$params$id),]
  shared_guess <- lapply(g$shared_params, function(x) runif(1, .3, 3.))
  
  res <- fitSharedParameters (shared_guess,
                       g$data,
                       forms,
                       g$params,
                       options,
                       g$size)
  abs(1 - unlist(g$shared_params) / unlist(res))
}

testFitDispersion <- function(n = 2, replicates = 2) {
  g <- generateTestData(n = n, replicates = replicates)
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    lower_boundary_size = 1 / 10,
    upper_boundary_size = 1e10,
    cores = 2
  )
  dispersion_guess <- runif(1, 1 / 10, 1e3)
  res <- fitDispersion(g$shared_params,
                       g$data,
                       forms,
                       g$params,
                       options,
                       dispersion_guess)
  abs(1 - g$size / res)
}

testFitModel <- function(n = 2, replicates = 2) {
  g <- generateTestData(n = n, replicates = replicates)
  d <- g$data
  d <- d[sample(length(d$id)), ]
  individual_params <- g$params[, which(names(g$params) != "id")]
  individual_params[, ] <- .1
  individual_params$id <- g$params$id
  options <- list(
    lower_boundary =c(1,1,0,0),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    lower_boundary_size = 1 / 10,
    upper_boundary_size = 1e10,
    cores = 2
  )
  guess <- guess_params(split(d, d$id), individual_params)
  shared_guess <- lapply(g$shared_params, function(x) runif(1, .3, 3.))
  fitResult <- fitModel(d,  forms, guess, shared_guess, options)
  p <- fitResult$individual_params
  g$params <- g$params[order(g$params$id),]
  p <- p[order(p$id),]
  errors <- (1 - p[, which(names(p) != "id"), drop = FALSE] /
                  g$params[, which(names(g$params) != "id"), drop = FALSE])
  list(
    individual_err = errors,
    shared_err = (1 - unlist(fitResult$shared_params) /
                       unlist(g$shared_params)),
    size_err = (1 - fitResult$size / g$size)
  )
}

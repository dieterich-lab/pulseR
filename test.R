#source("likelihoods.R")

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
  set.seed(259)
  forms <- getFormulas()
  conditions <- list()
  conditions$sample <- replicate(length(forms) * replicates,
                                 paste0(letters[sample(25, 10)], collapse = ""))
  conditions$condition <- names(forms)
  conditions <- as.data.frame(conditions, stringAsFactors = FALSE)
  genes <-
    replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- paste0("ENS00000", genes)
  p <- data.frame(
    mu_n = runif(n, 1e2, 50000),
    mu_h = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8),
    a_h  = runif(n, .05, .8)
  )
  rownames(p) <- genes
  alphas <- list(
    alpha_chase = 1,
    alpha_lab = 2,
    beta_chase = 1,
    beta_lab = 2
  )
  size <- 1e2
  d <- lapply(seq_along(rownames(p)), function(i) {
    generateTestDataSingle(forms, p[i,], alphas,
                           conditions = conditions$condition)
  })
  data <- do.call(rbind, d)
  rownames(data) <- genes
  colnames(data) <- conditions$sample
  rownames(conditions) <- conditions$sample
  list(
    count_data = data[order(rownames(data)), ],
    conditions = conditions[, -1, drop = FALSE],
    params = p[order(rownames(p)), ],
    size = size,
    shared_params = alphas
  )
}

getFormulas <- function() {
  MeanFormulas(
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
}

guess_params <- function(data, conditions) {
  totals <- data[, conditions$condition == "total_Norm", drop = FALSE]
  mean_expression <- apply(X = totals, FUN = mean, MARGIN = 1)
  guess <- data.frame(
    mu_n = mean_expression,
    mu_h = mean_expression,
    a_n = .1,
    a_h = .1
  )
  rownames(guess) <- rownames(data)
  guess
}

cookWorkEnvironment <- function(n, replicates) {
  formulas <- getFormulas()
  g <- generateTestData(n = n, replicates = replicates)
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    lower_boundary_size = 0,
    upper_boundary_size = 1e3,
    cores = 4
  )
  options$parscales <- mapply(max,
                           abs(options$upper_boundary),
                           abs(options$lower_boundary))
  pd <- PulseData(g$count_data, g$conditions, formulas)
  normalise(pd)
  g$count_data <- NULL
  g$conditions <- NULL
  list(pd = pd,
       options = options,
       params = g)
}

testIndividualGeneParams <- function(n = 2, replicates = 2) {
  wenv <- cookWorkEnvironment(n, replicates)
  pd <- wenv$pd
  g <- wenv$params
  guess <- guess_params(pd$count_data, pd$conditions)
  estimation <- fitIndividualParameters(
    old_params = guess,
    pulseData = pd,
    shared_params = g$shared_params,
    options = wenv$options,
    size = g$size
  )
  errors <- abs(1 - estimation[rownames(g$params), ] / g$params)
  errors
}

testSharedParams <- function(n = 2, replicates = 2) {
  wenv <- cookWorkEnvironment(n, replicates)
  pd <- wenv$pd
  g <- wenv$params
  shared_guess <-
    lapply(g$shared_params, function(x)
      runif(1, .3, 3.))
  res <- fitSharedParameters (
    old_shared_params = shared_guess,
    pulseData = pd,
    individual_params = g$params,
    options = wenv$options,
    size = g$size
  )
  abs(1 - unlist(g$shared_params) / unlist(res))
}

testFitDispersion <- function(n = 2, replicates = 2) {
  wenv <- cookWorkEnvironment(n, replicates)
  pd <- wenv$pd
  g <- wenv$params
  dispersion_guess <- runif(1, 1 / 10, 1e3)
  res <- fitDispersion(
    shared_params = g$shared_params,
    pulseData = pd,
    individual_params = g$params,
    options = wenv$options,
    size = dispersion_guess
  )
  abs(1 - g$size / res)
}

testFitModel <- function(n = 2, replicates = 2) {
  wenv <- cookWorkEnvironment(n, replicates)
  pd <- wenv$pd
  g <- wenv$params
  guess <- guess_params(pd$count_data, pd$conditions)
  shared_guess <-
    lapply(g$shared_params, function(x)
      runif(1, .3, 3.))
  fitResult <- fitModel(
    pulseData = pd,
    params        = guess,
    shared_params = shared_guess,
    options       = wenv$options
  )
  p <- fitResult$individual_params
  errors <- (1 - p / g$params)
  list(
    individual_err = errors,
    shared_err = (
      1 - unlist(fitResult$shared_params) /
        unlist(g$shared_params)
    ),
    size_err = (1 - fitResult$size / g$size)
  )
}

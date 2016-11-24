source("likelihoods.R")

generateTestDataSingle <- function(forms,
                                   individual_params,
                                   shared_params,
                                   n = 1,
                                   size = 100) {
  means <- sapply(forms, eval,
                  c(as.list(individual_params), as.list(shared_params)))
  counts <- rnbinom(n = length(means) * n,
                    mu = means,
                    size = size)
  data.frame(
    condition = rep(names(forms), n),
    count = counts,
    norm_factor = 1
  )
}

generateTestData <- function(n, replicates) {
  genes <- letters[1:n]
  p <- data.frame(
    id = genes,
    mu_n = runif(n, 1e3, 5000),
    mu_h = runif(n, 1e3, 5000),
    a_n  = runif(n, .05, .8),
    a_h  = runif(n, .05, .8)
  )
  alphas <- list(
    alpha_chase = 1,
    alpha_lab = 1,
    beta_chase = 1,
    beta_lab = 1
  )
  size <- 1e2
  d <- lapply (seq_along(p$id), function(i) {
    cbind(id = p$id[i],
          generateTestDataSingle(forms, p[i, -1], alphas, n = replicates))
  })
  data <- do.call(rbind, d)
  list(data = data,
       params = p,
       shared_params = alphas)
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

testIndividualGeneParams <- function(n = 2, replicates = 2) {
  g <- generateTestData(n = n, r = replicates)
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4),
    cores = 2
  )
  data <- split(g$data, g$data$id)
  mean_expression <-
    unlist(lapply(data, function(x)
      mean(x$count[x$condition == "total_Norm"])))
  guess <- data.frame(
    id = g$params$id,
    mu_n = mean_expression,
    mu_h = mean_expression,
    a_n = .1,
    a_h = .1
  )
  estimation <-
    fitIndividualParameters(guess, data, forms, g$shared_params, options, size)
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
  f <- ll_shared_params(
    count_data = g$data,
    forms = forms,
    individual_params =  g$params,
    shared_param_names = names(g$shared_params),
    size = 1e2
  )
  res <- optim(
    par = runif(4, .3, 3),
    fn = f,
    method = "L-BFGS-B",
    lower = rep(1e-9, length(g$shared_params)),
    upper = c(15, 15, 4, 4)
  )
  abs(1 - unlist(g$shared_params) / res$par)
}

testFitModel <- function() {
  # Make data
  d <- list()
  p <- data.frame(
    mu_n = c(100, 1000),
    mu_h = c(50, 500),
    a_n = c(.5, .8),
    a_h = c(.8, .5)
  )
  param_names <- names(p)
  alphas <-
    list(
      alpha_chase = 2,
      alpha_lab = 1.5,
      beta_chase = 1,
      beta_lab = .8
    )
  p$id <- c("b", "a")
  for (i in seq_along(p$id)) {
    d[[i]] <- cbind(id = p$id[i],
                    generateTestData(forms, p[i, param_names], alphas, n = 2))
  }
  d <- do.call(rbind, d)
  d$norm_factor <- 1
  d <- d[sample(length(d$id)),]
  # Fit
  options <- list(
    lower_boundary = rep(1e-9, 4),
    upper_boundary = c(1e5, 1e5, 1, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 4),
    upper_boundary_shared = rep(5, 4)
  )
  individual_params <- as.data.frame(p[, -5])
  individual_params[,] <- .1
  individual_params$id <- p$id
  fitResult <-
    fitModel (d,  forms, individual_params, alphas, options)
  ip <- rbind(estimated = fitResult$individual_params,
              correct = p)
  ip <- split(ip, ip$id)
  sp <- rbind(estimated = fitResult$shared_params, correct = alphas)
  print(ip)
  print(sp)
}

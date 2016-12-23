source("likelihoods.R")
source("fit.R")
source("pulseData.R")

# Create a test data set according to formulas form
# parameters par${individual_params, shared_params, fraction_factors}
# conditions${condition, fraction}
# Returns a matrix with rows ordered as par$individual_params and
# columns ordered as conditions row
generateTestDataFrom <- function(forms, par, conditions) {
  counts <- list()
  for(i in seq_along(par$individual_params[,1])){
    means <- sapply(forms, eval, 
      c(as.list(par$individual_params[i,]),
        as.list(par$shared_params)))
    # normalise
    if(!is.null(conditions$fraction_indexes)){
      fraction_indexes <- match(as.integer(conditions$fraction))
      means <- means * c(1, par$fraction_factors)[fraction_indexes]
    }
    indexes <- match(conditions$condition, names(forms))
    counts[[i]] <- rnbinom(
      n    = length(conditions$condition),
      mu   = means[indexes],
      size = par$size)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- rownames(par$individual_params)
  colnames(counts) <- rownames(conditions)
  counts
}

conditionsFromFormulas <- function(forms, replicates) {
  conditions <-rep(names(forms), replicates)
  names(conditions) <- replicate(length(forms) * replicates,
    paste0("sample_", paste0(letters[sample(25, 10)], collapse = "")))
  conditions
}

generateTestData <- function(n,
                             replicates,
                             forms,
                             conditions){
  genes <- replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- paste0("ENS00000", genes)
  par <- list()
  par$individual_params <- data.frame(
    mu_n = runif(n, 1e2, 50000),
    mu_h = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8),
    a_h  = runif(n, .05, .8)
  )
  rownames(par$individual_params) <- genes
  par$individual_params <- par$individual_params[order(rownames(par$individual_params)),]
  par$shared_params <- list(
    alpha_chase = 1,
    alpha_lab = 2,
    beta_chase = 1,
    beta_lab = 2
  )
  par$size <- 1e2
  d <- generateTestDataFrom(forms, par, conditions)
  d <- d[order(rownames(d)),]
  list( count_data=d, par=par, conditions = conditions)
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

getFormulasWithHyperParams <- function() {
  MeanFormulas(
    total = mu_n,
    flow_lab      = alpha_lab * (mu_n * a_n^time),
    biotin_lab    = beta_lab * mu_n * (1 - a_n^time)
  )
}

guess_params <- function(wenv) {
  guess <-(apply(wenv$par$individual_params, 2,median))
  guess <- matrix(guess, nrow=dim(wenv$par$individual_params)[1],
                  ncol=length(guess), byrow=TRUE)
  rownames(guess) <- rownames(wenv$par$individual_params)
  colnames(guess) <- colnames(wenv$par$individual_params)
  guess + 1e-10
}

cookWorkEnvironment <- function(n,
                                replicates,
                                formulas=getFormulas(),
                                conditions) {
  set.seed(259)
  conditions <- list()
  conditions$condition <- conditionsFromFormulas(forms = formulas,
                                       replicates = replicates)
  conditions <- as.data.frame(conditions)
  g <- generateTestData(n = n,
                        replicates = replicates,
                        forms = formulas, 
                        conditions = conditions)
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
       par = g$par)
}

testIndividualGeneParams <- function(wenv, thres=.05) {
  pd <- wenv$pd
  true_params <- wenv$par$individual_params
  par <- wenv$par
  par$individual_params <- guess_params(wenv)
  estimation <- fitIndividualParameters(
    pulseData = pd,
    par = wenv$par,
    options = wenv$options
  )
  estimation <- estimation[rownames(wenv$par$individual_params), ]
  errors <- abs(1 - estimation / true_params)
  stopifnot(max(errors) < thres)
  errors
}

testSharedParams <- function(wenv, thres=0.05) {
  pd <- wenv$pd
  true_params <- wenv$par$shared_params
  shared_guess <- lapply(wenv$par$shared_params, function(x) runif(1, .3, 3.))
  wenv$par$shared_params <- shared_guess
  res <- fitSharedParameters (pd, wenv$par, wenv$options)
  errors <- abs(1 - unlist(true_params) / unlist(res))
  stopifnot(max(errors) < thres)
  errors
}

testFitDispersion <- function(wenv, thres=0.05) {
  pd <- wenv$pd
  dispersion_guess <- runif(1, 1 / 10, 1e3)
  res <- fitDispersion(pd, wenv$par,  wenv$options)
  errors <- abs(1 - wenv$par$size / res)
  stopifnot(max(errors) < thres)
  errors
}

testFitModel <- function(wenv, thres=0.05) {
  pd <- wenv$pd
  par <- wenv$par
  par$size <- 10
  par$individual_params <- guess_params(wenv)
  shared_guess <- lapply(par$shared_params, function(x) runif(1, .3, 3.))
  par$shared_params <- shared_guess
  fitResult <- fitModel(pd, par, wenv$options)
  p <- fitResult$par$individual_params
  p <- p[rownames(par$individual_params), ]
  errors <- (1 - p /wenv$par$individual_params)
  res <- list(
    individual_err = errors,
    shared_err = ( 1 - unlist(fitResult$par$shared_params) / unlist(wenv$par$shared_params)),
    size_err = (1 - fitResult$par$size / wenv$par$size)
  )
  stopifnot(max(unlist(res)) < thres)
  res
}


testAll <- function(n=50, replicates=50){
  wenv <- cookWorkEnvironment(n, replicates,  getFormulas()) 
  testIndividualGeneParams(wenv) 
  print("individual OK")
  testSharedParams(wenv) 
  print("shared OK")
  testFitDispersion(wenv)
  print("dispersion OK")
  #testFitModel(wenv)
  return("OK")
}

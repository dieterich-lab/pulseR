## Tests for data with several time points
source("test.R")
cookWorkEnvironmentWithTime <- function(n,
                                replicates,
                                time_n = 3,
                                formulas=getFormulasWithHyperParams()) {
  set.seed(259)
  conditions <-data.frame(
    condition = conditionsFromFormulas(forms = formulas, 
                                       replicates = replicates * time_n))
  conditions$time <- rep(1:time_n, each=length(formulas))
  conditions$dummy <- "i_will_fail_your_code"
  t <- addKnownShared(formulas, conditions)
  conditions <- t$conditions
  formulas <- t$formulas
  g <- generateTestDataWithTime(n = n,
                        replicates = replicates,
                        forms = formulas,
                        conditions = conditions)
  options <- list(
    lower_boundary = rep(1e-9, 2),
    upper_boundary = c(1e5, 1) - 1e-1,
    lower_boundary_shared = rep(1e-9, 2),
    upper_boundary_shared = rep(5, 2),
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

generateTestDataWithTime <- function(n,
                             replicates,
                             forms,
                             conditions){
  genes <- replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- paste0("ENS00000", genes)
  p <- data.frame(
    mu_n = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8)
  )
  rownames(p) <- genes
  alphas <- list(
    alpha_lab = 1,
    beta_lab = 2
  )
  size <- 1e2
  d <- lapply(seq_along(rownames(p)), function(i) {
    generateTestDataSingle(forms, p[i,], alphas,
                           conditions = conditions)
  })
  data <- do.call(rbind, d)
  rownames(data) <- genes
  colnames(data) <- names(conditions)
  list(
    count_data = data[order(rownames(data)), ],
    conditions = conditions,
    params = p[order(rownames(p)), ],
    size = size,
    shared_params = alphas
  )
}

testTimeData <- function(n=10, replicates=10){
  wenv <- cookWorkEnvironmentWithTime (
    n, replicates, time_n = 3, getFormulasWithHyperParams()) 
  testIndividualGeneParams(n, replicates, wenv) 
  testSharedParams(n, replicates, wenv) 
  testFitModel(n, replicates, wenv)
  return("OK")
}

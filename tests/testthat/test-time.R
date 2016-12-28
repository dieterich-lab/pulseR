## Tests for data with several time points

getFormulasWithHyperParams <- function() {
  MeanFormulas(
    total = mu_n,
    flow_lab      = (mu_n * a_n^time),
    biotin_lab    =  mu_n * (1 - a_n^time)
  )
}

cookWorkEnvironmentWithTime <- function(n,
                                replicates,
                                time_n = 3,
                                formulas=getFormulasWithHyperParams()) {
  set.seed(259)
  conditions <-data.frame(
    condition = conditionsFromFormulas(forms = formulas, 
                                       replicates = replicates * time_n))
  conditions$fraction <- conditions$condition # add scale info
  conditions$fraction <- relevel(conditions$fraction, "total")
  conditions$time <- rep(1:time_n, each=length(formulas))
  conditions$dummy <- "i_will_fail_your_code"
  t <- addKnownShared(formulas, conditions)
  g <- generateTestDataWithTime(n  = n,
                        replicates = replicates,
                        formulas   = formulas,
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
  pd <- PulseData(
    count_data = g$count_data,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
  normalise(pd) 
  g$count_data <- NULL
  g$conditions <- NULL
  list(pd = pd,
       options = options,
       par = g$par)
}

generateTestDataWithTime <- function(n, replicates, formulas, conditions){
  genes <- replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- sort(paste0("ENS00000", genes))
  p <- data.frame(
    mu_n = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8)
  )
  rownames(p) <- genes
  par <- list()
  par$individual_params <- p
  par$fraction_factors <- c(1,1)
  par$size <- 1e2
  t <- addKnownShared(formulas, conditions)
  conditions$condition <- t$conditions
  formulas <- t$formulas
  data <- generateTestDataFrom(formulas, par, conditions)
  list(
    count_data = data,
    par=par,
    conditions = conditions
  )
}

testTimeData <- function(n=10, replicates=10, thres=.05){
  wenv <- cookWorkEnvironmentWithTime (
    n, replicates, time_n = 3, getFormulasWithHyperParams()) 
  testIndividualGeneParams(wenv, thres) 
  testFitModel(wenv, thres)
  return("OK")
}

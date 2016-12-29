
# Create a test data set according to formulas form
# parameters par${individual_params, shared_params, fraction_factors}
# conditions${condition, fraction}
# Returns a matrix with rows ordered as par$individual_params and
# columns ordered as conditions row
conditionsFromFormulas <- function(forms, replicates) {
  conditions <-rep(names(forms), replicates)
  names(conditions) <- sort(replicate(length(forms) * replicates,
    paste0("sample_", paste0(letters[sample(25, 10)], collapse = ""))))
  conditions
}

# Create a test data set from random parameters
# no fraction or time information is used
generateTestData <- function(n, replicates, forms, conditions){
  genes <- replicate(n, paste0(letters[sample(25, 10)], collapse = ""))
  genes <- sort(paste0("ENS00000", genes))
  par <- list()
  par$individual_params <- data.frame(
    mu_n = runif(n, 1e2, 50000),
    mu_h = runif(n, 1e2, 50000),
    a_n  = runif(n, .05, .8),
    a_h  = runif(n, .05, .8)
  )
  rownames(par$individual_params) <-genes
  par$shared_params <- list(
    alpha_chase = 1,
    alpha_lab = 2,
    beta_chase = 1,
    beta_lab = 2
  )
  par$size <- 1e2
  d <- generateTestDataFrom(forms, par, conditions)
  list( count_data=d, par=par, conditions = conditions)
}

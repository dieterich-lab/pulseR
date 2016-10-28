source("likelihoods.R") 
library(testthat)

generateTestData <- function(forms, individual_params, shared_params, n=1){
   means <- sapply(forms, eval, c(as.list(individual_params), as.list(shared_params)))
   counts <- rpois(n=length(means)*n, lambda=means)
   data.frame(condition=rep(names(forms), n), count=counts)
}

forms <- list(
     total_Hypox        = quote(mu_n),
     total_Norm         = quote(mu_n*a_h + mu_h),
     flow_lab_Norm      = quote(alpha_lab*(mu_n*a_n)),
     biotin_lab_Norm    = quote(beta_lab*mu_n*(1-a_n)),
     flow_lab_Hypox     = quote(alpha_lab*(mu_n*a_h)),
     biotin_lab_Hypox   = quote(beta_lab*mu_h),
     flow_chase_Norm    = quote(alpha_chase*mu_n*(1-(1-a_n )*a_n )),
     biotin_chase_Norm  = quote(beta_chase*mu_n*(1-a_n )*a_n),
     flow_chase_Hypox   = quote(alpha_chase*(mu_h + mu_n*a_n * a_h)),
     biotin_chase_Hypox = quote(beta_chase*(mu_n*(1-a_n)*a_h)))
     
p <- c(mu_n=1000, mu_h=500, a_n=.1, a_h=.2)
alphas <- list(alpha_chase=1, alpha_lab=1, beta_chase=1,beta_lab=1)
d <- generateTestData(forms, p,alphas, n=2)
d$deseq_factor <- 1

likelihood <- ll_gene(d$count, d$deseq_factor, d$condition,
    forms, param_names=names(p), alphas )

res <-optim(p, likelihood, method="L-BFGS-B", 
    lower=rep(1e-9, length(params)), upper=c(1e5,1e5, 1,1))
expect_equal(res$par,p)
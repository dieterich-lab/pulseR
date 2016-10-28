source("likelihoods.R") 

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
     
testIndividualGeneParams <- function(){
    p <- c(mu_n=1000, mu_h=500, a_n=.1, a_h=.2)
    alphas <- list(alpha_chase=1, alpha_lab=1, beta_chase=1,beta_lab=1)
    d <- generateTestData(forms, p,alphas, n=2)
    d$deseq_factor <- 1

    likelihood <- ll_gene(d$count, d$deseq_factor, d$condition,
        forms, param_names=names(p), alphas )

    res <-optim(p, likelihood, method="L-BFGS-B", 
        lower=rep(1e-9, length(p)), upper=c(1e5,1e5, 1,1))
    results <- cbind(correct=p,estimated=res$par)
    rownames(results)
    results
}

testSharedParams <- function(){
    d <- list()
    p <- data.frame(mu_n=c(100,1000), 
                    mu_h=c(50,500), 
                     a_n=c(.5, .8),
                     a_h=c(.8, .5))
    alphas <- list(alpha_chase=2, alpha_lab=1.5, beta_chase=1,beta_lab=.8)
    p$id <- c("a", "b")
    for(i in seq_along(p$id)){
        d[[i]] <- cbind(id=p$id[i],generateTestData(forms, p[i,1:4],alphas, n=2))
    }
    d <- do.call(rbind,d)
    d$deseq_factor <- 1
    f <- ll_shared_params (d, forms,p, names(alphas))
    res <-optim(rep(1,4), f, method="L-BFGS-B", 
        lower=rep(1e-9, length(alphas)), upper=c(15,15, 4,4))
    res <- cbind(correct=unlist(alphas), estimated=res$par)
    rownames(res) <- names(alphas)
    res
}




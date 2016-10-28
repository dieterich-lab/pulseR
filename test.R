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
    d$norm_factor<- 1

    likelihood <- ll_gene(d, forms, param_names=names(p), alphas )

    res <-optim(p, likelihood, method="L-BFGS-B", 
        lower=rep(1e-9, length(p)), upper=c(1e5,1e5, 1,1))
    results <- cbind(correct=p,estimated=res$par)
    rownames(results)
    results
}

# > testIndividualGeneParams()
#      correct    estimated
# mu_n   1e+03 998.12154366
# mu_h   5e+02 491.52358247
# a_n    1e-01   0.09599496
# a_h    2e-01   0.19939628

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
    d$norm_factor<- 1
    f <- ll_shared_params (d, forms,p, names(alphas))
    res <-optim(rep(1,4), f, method="L-BFGS-B", 
        lower=rep(1e-9, length(alphas)), upper=c(15,15, 4,4))
    res <- cbind(correct=unlist(alphas), estimated=res$par)
    rownames(res) <- names(alphas)
    res
}

# > testSharedParams()
#             correct estimated
# alpha_chase     2.0 1.9863519
# alpha_lab       1.5 1.5321680
# beta_chase      1.0 0.9276926
# beta_lab        0.8 0.7643754



testFitModel <- function(){
    # Make data
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
    d$norm_factor<- 1
    # Fit
    lower_boundary <- rep(1e-9,4)
    upper_boundary <- c(1e5,1e5,1,1)
    lower_boundary_shared <- rep(1e-9,4)
    upper_boundary_shared <- rep(5)
    individual_params <- as.data.frame(p[,-5])
    individual_params[,] <- 1
    individual_params$id <- p$id
    fitResult <- fitModel (d,  forms, individual_params,
                alphas,
                lower_boundary, 
                upper_boundary,
                lower_boundary_shared,
                upper_boundary_shared)
    ip <- rbind(estimated=fitResult$individual_params, correct=p)
    ip <- split(ip, ip$id)
    sp <- rbind(estimated=fitResult$shared_params, correct=alphas)
    print(ip)
    print(sp)
}

# > testFitModel()
# $a
#                 mu_n     mu_h       a_n       a_h id
# estimated.a 101.2886 56.17576 0.4852201 0.7556706  a
# correct.1   100.0000 50.00000 0.5000000 0.8000000  a
# 
# $b
#                 mu_n     mu_h       a_n       a_h id
# estimated.b 1009.863 497.4247 0.7981516 0.4971237  b
# correct.2   1000.000 500.0000 0.8000000 0.5000000  b
# 
#           alpha_chase alpha_lab beta_chase beta_lab 
# estimated 1.992333    1.503761  0.9536717  0.8084433
# correct   2           1.5       1          0.8  
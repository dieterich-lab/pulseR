#from pryr package
 substitute_q <- function (x, env) 
{
    stopifnot(is.language(x))
    env <- to_env(env)
    call <- substitute(substitute(x, env), list(x = x))
    eval(call)
}

# from pryr package
to_env <-  function (x, quiet = FALSE) 
{
    if (is.environment(x)) {
        x
    }
    else if (is.list(x)) {
        list2env(x)
    }
    else if (is.function(x)) {
        environment(x)
    }
    else if (length(x) == 1 && is.character(x)) {
        if (!quiet) 
            message("Using environment ", x)
        as.environment(x)
    }
    else if (length(x) == 1 && is.numeric(x) && x > 0) {
        if (!quiet) 
            message("Using environment ", search()[x])
        as.environment(x)
    }
    else {
        stop("Input can not be coerced to an environment", call. = FALSE)
    }
}

 
makeVector <- function(forms){
    string <- paste("c(", 
        paste(mapply(paste,names(forms), as.character(forms), sep="="), 
            collapse=","),
            ")") 
    parse(text=string)
}


ll_gene <- function(count_data, forms, param_names, shared_params=NULL){
    conditions <- count_data$condition
    counts <- count_data$count
    norm_factors <- count_data$norm_factor
    mean_indexes <- sapply(conditions, match, names(forms))
    if(!is.null(shared_params)&&!anyNA(forms))
        forms <- lapply(forms, substitute_q, shared_params)
    means_vector<- makeVector(forms)
    funquote <- substitute(
        function(params, output=FALSE){
            names(params) <- param_names
            mus <- eval(means_vector, as.list(params))
            lambdas <- mus[mean_indexes]
            if(output){
                return(data.frame(lambdas=lambdas*norm_factors, count=counts))
            }
            -sum(dpois(counts, (lambdas+1e-10)* norm_factors, log=TRUE))
        }, parent.frame()
    )
    eval(funquote)
}


ll_shared_params <- function(count_data,forms,individual_params,
                             shared_param_names){
    means_for_genes <- list()
    individual_params$id <- as.factor(individual_params$id)
    individual_params <- individual_params[order(individual_params$id),]
    forms_vector <- makeVector(forms)[[1]]
    condition_indexes <- sapply(count_data$condition, match, names(forms))
    means_indexes <- cbind(as.numeric(count_data$id), condition_indexes)
    for(i in 1:2){
        means_for_genes[[i]] <- substitute_q(
            forms_vector, env=individual_params[i,])
    }
    function(shared_params){
        names(shared_params ) <- shared_param_names
        means_matrix <- matrix(0,ncol=length(forms), nrow=2)
        for(i in 1:2){
            means_matrix[i,] <- eval(means_for_genes[[i]],
                                     as.list(shared_params))
        }
        lambdas <- means_matrix[means_indexes]
        -sum(dpois(count_data$count, (lambdas+1e-10), log=TRUE))
    }
}


fitModel <- function(count_data,  formulas, individual_params,
                     shared_params, lower_boundary, upper_boundary){
    # Fit params for every genes individually
    splitted_data <- split(count_data, count_data$id)
    param_names <- names(individual_params)
    param_names <- param_names[-which(param_names=="id")]
    new_params <- list()
    for(gene in names(splitted_data)){
        objective <- ll_gene(splitted_data[[gene]], forms, 
                             param_names, shared_params)
        new_params [[gene]] <- optim(rep(1,4), objective, method="L-BFGS-B",
            lower=lower_boundary, upper=upper_boundary)$par
    }
    individual_params <- as.data.frame(do.call(rbind, new_params))
    names(individual_params) <- param_names
    individual_params$id <- rownames(individual_params)
    individual_params
    # Fit shared params
    
}

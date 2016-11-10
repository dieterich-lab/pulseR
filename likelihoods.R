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

MeanFormulas <- function(...){
    eval(substitute(alist(...)))
}
 
makeVector <- function(forms){
    string <- paste("c(", 
        paste(mapply(paste,names(forms), as.character(forms), sep="="), 
            collapse=","),
            ")") 
    parse(text=string)
}

contaminate <- function(formulas, target_condition, 
                        contaminant_condition, coef_name){
    target_condition <- deparse(substitute(target_condition))
    contaminant_condition <- deparse(substitute(contaminant_condition))
    f1 <- deparse(formulas[[target_condition]])
    f2 <- deparse(formulas[[contaminant_condition]])
    e <- paste("(1-",coef_name,")*(",f1,")+",coef_name,"*(",f2,")")
    parse(text=e)[[1]]
}


ll_gene <- function(count_data, forms, param_names, shared_params=NULL){
    conditions <- count_data$condition
    counts <- count_data$count
    norm_factors <- count_data$norm_factor
    stopifnot(!is.null(conditions),!is.null(norm_factors),!is.null(counts))
    mean_indexes <- sapply(conditions, match, names(forms))
    if(!is.null(shared_params)&&!anyNA(forms))
        forms <- lapply(forms, substitute_q, shared_params)
    means_vector<- makeVector(forms)
    funquote <- substitute(
        function(params, output=FALSE){
            names(params) <- param_names
            mus <- eval(means_vector, as.list(params))
            lambdas <- mus[mean_indexes]
            -sum(dpois(counts, (lambdas+1e-10)* norm_factors, log=TRUE))
        }, parent.frame()
    )
    eval(funquote)
}

getMeansEstimatingFunction <- function(count_data, individual_params, 
                              forms, shared_param_names){
    gene_number <- length(individual_params$id)
    individual_params$id <- as.factor(individual_params$id)
    individual_params <- individual_params[order(individual_params$id),]
    forms_vector <- makeVector(forms)[[1]]
    condition_indexes <- sapply(count_data$condition, match, names(forms))
    means_indexes <- cbind(as.numeric(count_data$id), condition_indexes)
    means_for_genes <- list()
    for(i in 1:gene_number){
        means_for_genes[[i]] <- substitute_q(
            forms_vector, env=individual_params[i,])
    }
    means_matrix <- matrix(0,ncol=length(forms), nrow=gene_number)
    function(shared_params){
        names(shared_params ) <- shared_param_names
        for(i in 1:gene_number){
            means_matrix[i,] <- eval(means_for_genes[[i]],
                                     as.list(shared_params))
        }
        means_matrix[means_indexes]
    }
}

ll_shared_params <- function(count_data,forms,individual_params,
                             shared_param_names){
    estimateMeans <- getMeansEstimatingFunction(count_data, individual_params, 
                                        forms, shared_param_names)
    function(shared_params){
        lambdas <- estimateMeans(shared_params)
        #-sum(dpois(count_data$count, (lambdas+1e-10), log=TRUE))
        -median(dpois(count_data$count, (lambdas+1e-10), log=TRUE))
    }
}

predict.expression <- function(count_data, model, forms){
    estimateMeans <- getMeansEstimatingFunction(count_data,
    model$individual_params, forms, names(model$shared_params))
    lambdas <- estimateMeans(model$shared_params)
    list(prediction=lambdas,
         logL=dpois(count_data$count, (lambdas+1e-10), log=TRUE))
}


fitIndividualParameters <- function(old_params, splitted_data, formulas,
                                    shared_params, options){
    verbose <- options$verbose
    if(verbose=="verbose"){
        ngene <- dim(old_params)[1]
    }
    ids <- old_params$id
    param_names <- setdiff(names(old_params), "id")
    old_params <- split(old_params[,param_names],
                        old_params$id)
    new_params <- list()
    for(gene_index in seq_along(old_params)){
        objective <- ll_gene(splitted_data[[gene_index]], formulas, 
                            param_names, shared_params)
        new_params [[gene_index]] <- optim(
            unlist(old_params[[gene_index]]), 
            objective, 
            method="L-BFGS-B", 
            lower=options$lower_boundary, 
            upper=options$upper_boundary)$par
        if(verbose=="verbose"){
            cat(rep("",100),"\r")
            cat(gene_index, " of ", ngene, " are analysed\r")
        }
    }
    new_params <- as.data.frame(do.call(rbind, new_params))
    names(new_params) <- param_names
    new_params$id <- ids
    new_params
}

fitSharedParameters <- function(old_shared_params, count_data, formulas,
                                individual_params, options){
    shared_objective <- ll_shared_params(count_data, formulas, 
        individual_params, names(old_shared_params))
    shared_params <- optim(
        unlist(old_shared_params),
        shared_objective, 
        method="L-BFGS-B", 
        lower=options$lower_boundary_shared, 
        upper=options$upper_boundary_shared)$par
    names(shared_params) <- names(old_shared_params)
    as.list(shared_params)
}

# options is a list with records
# - individual_rel_err
# - shared_rel_tol
fitModel <- function(count_data,  formulas, individual_params,
                     shared_params=NULL, options=list()){
    splitted_data <- split(count_data, count_data$id)
    param_names <- setdiff(names(individual_params), "id")
    shared_param_names <- names(shared_params)
    opts <- list(
        individual_rel_tol=rep(1e-2, length(param_names)),
        shared_rel_tol=rep(1e-2, length(shared_param_names)),
        verbose="silent")
    individual_rel_err <- 10*opts$individual_rel_tol
    shared_rel_err <- 10* opts$shared_rel_tol
    if(is.null(shared_params)) shared_rel_err <- 0
    while(any(individual_rel_err > opts$individual_rel_tol) || 
        any(shared_rel_err > opts$shared_rel_tol)){
        shared_params <- as.list(shared_params)
        opts[names(options)] <- options
        # Fit params for every genes individually
        old_params <- individual_params
        individual_params <- fitIndividualParameters(
            old_params, splitted_data, formulas, shared_params, opts)
        individual_rel_err <- max(
            abs(1 - individual_params[,param_names] / old_params[,param_names]))
        # Fit shared params
        if(!is.null(shared_params)){
            old_shared_params <- shared_params
            shared_params <- fitSharedParameters(shared_params, count_data,
                formulas, individual_params, opts)
            shared_rel_err <- 
                1 - unlist(shared_params) / unlist(old_shared_params)
        }
    }
    list(individual_params=individual_params, shared_params=shared_params)
}

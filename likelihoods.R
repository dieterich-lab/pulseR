 
makeVector <- function(forms){
    string <- paste("c(", 
        paste(mapply(paste,names(forms), as.character(forms), sep="="), 
            collapse=","),
            ")") 
    parse(text=string)
}


ll_gene <- function(counts, norm_factors, conditions, 
                    forms, param_names, param_list=NULL){
    mean_indexes <- sapply(conditions, match, names(forms))
    means_vector<- makeVector(forms)
    funquote <- substitute(
        function(params, output=FALSE){
            names(params) <- param_names
            env <- c(as.list(params),param_list )
            mus <- eval(means_vector, env)
            lambdas <- mus[mean_indexes]
            if(output){return(data.frame(lambdas=lambdas*norm_factors, count=counts))}
            -sum(dpois(counts, (lambdas+1e-10)* norm_factors, log=TRUE))
        } 
    )
    eval(funquote)
}


forms <- alist(
  lab = a,
  unlab = b,
  total = b + a)
  
counts <- data.frame(
  lab=1:10, 
  unlab=runif(10)*10,
  total = abs(rnorm(10)) +1:10)
  
par <- data.frame(a=1:10, b=10:1)

form_indexes <- list(
  c(1,2),
  c(1,2),
  3
)

norm_factors <- list(
  c(5,6),
  c(2,7),
  1
)

llgene <- function(params){
  
}

evaled <- lapply(forms, eval, par)

sample_means <- function(evaled_forms, form_indexes, norm_factors){
  mus <- mapply(
    function(i, n){
      m <- do.call(cbind,evaled_forms[i])
      m %*% n
    },
    form_indexes,
    norm_factors, SIMPLIFY=FALSE)
  do.call(cbind, mus)
}

ll <- function(mus, counts, size){
  sum(dnbinom(counts, mus, log=TRUE, size=size))
}

.mu <- function(evaled_form_list, norm_fac){
  
}


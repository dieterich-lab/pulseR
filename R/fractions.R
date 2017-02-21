
forms <- alist(
  lab = a,
  unlab = b,
  total = b + a)
  
nGenes <- 7
## params are just env
par <- list(size = 1e2)
par$alpha <-  1
par <-  c(par, list(
  a = (1:nGenes) * 1e5, b = rep(1e6, nGenes)))
par$size <- 100000


counts <- as.matrix(data.frame(
  lab = rnbinom(nGenes, size = par$size, mu = par$a),
  unlab = rnbinom(nGenes, size = par$size, mu = par$b),
  total = rnbinom(nGenes, size = par$size, mu = par$a + par$b)
))
  

formIndexes <- list(
  c(1,2),
  c(1,2),
  3
)

normFactors <- list(
  c(1,0),
  c(0,1),
  1
)

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

ll <- function(counts, par, namesToOptimise) {
  p <- par[namesToOptimise]
  par[namesToOptimise] <- NULL
  function(x){
    par[namesToOptimise] <- relist(x, p)
    evaledForms <- lapply(forms, eval, envir = par) 
    means <- sample_means(evaledForms, formIndexes, normFactors)
    -sum(dnbinom(counts, mu = means, size = par$size, log = TRUE))
  }
}

f <- ll(counts, par, c("a", "b"))

paramClasses <- list(
  part = c("a", "b"),
  shared = c("alpha"))

# extend boundaries to param length
.b <- function(b, par) {
  for (p in names(b)) {
    if (length(b[[p]]) == 1)
      b[[p]] <- rep(b[[p]], length(par[[p]]))
  }
  b
}

fitParams <- function(counts, par, namesToOptimise, opts){
 objective <- ll(counts, par, namesToOptimise)
 lb <- unlist(opts$lb[namesToOptimise])
 ub <- unlist(opts$ub[namesToOptimise])
 x <- unlist(par[namesToOptimise])
 x <- optim(
   x,
   objective,
   method = "L-BFGS-B",
   control = list(parscale = x),
   lower = lb,
   upper = ub
 )$par
 relist(x, par[namesToOptimise])
}

opts <- list()
opts$lb <- list(a=rep(.1, b=.1))
opts$lb <- .b(opts$lb, par)
opts$ub <- list(a=1e7, b=2e7)
opts$ub <- .b(opts$ub, par)
                
fitParams(counts, par, c("a", "b"), opts)


forms <- alist(
  lab = a*time,
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
  s1 = rnbinom(nGenes, size = par$size, mu = par$a),
  s2 = rnbinom(nGenes, size = par$size, mu = par$b),
  s3 = rnbinom(nGenes, size = par$size, mu = par$a + par$b)
))
  
conditions <- data.frame(
  fraction = c("lab", "unlab", "total"),
  time = c(2,2,1))
pulseR:::addKnownShared(forms, conditions)
formulaIndexes <- list(
  lab = c("lab", "unlab"),
  unlab = c("lab", "unlab"),
  total = c('total'))



pd$formIndexes <- lapply(formulaIndexes, match, names(forms))  
  
pd <- list(counts = counts)
pd$formulas <- forms

pd$formIndexes <- list(
  c(1,2),
  c(1,2),
  3
)

pd$normFactors <- list(
  c(1,0),
  c(0,1),
  1
)

evaled <- lapply(forms, eval, par)


opts <- list()
opts$lb <- list(a=.1, b=.1)
opts$lb <- pulseR:::.b(opts$lb, par)
opts$ub <- list(a=1e7, b=2e7)
opts$ub <- pulseR:::.b(opts$ub, par)
                
fitParams(pd, par, c("a", "b"), opts)

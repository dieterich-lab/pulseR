# get index of the parameter on the `path` in the unlist(params)
.getIndex <- function(path, params) {
  x <- unlist(params)
  x[] <- seq_along(x)
  index <- .getElement(relist(x, params), path)
  if (length(index) > 1) 
    stop("Not complete path to the element (length of selection > 1)")
  index
}
  
# assign value in the list of parameters `params` given the path as a list
# `path` (e.g. list("mu", 2) will assign to params$mu[2])
.assignElement <- function(params, path, value) {
  if (length(path) == 1) {
    params[[unlist(path)]] <- value 
    params
  } else {
    params[[unlist(path[1])]] <-
      Recall(params[[unlist(path[1])]], path[-1], value)
    params 
  }
}

# get the element value located on the `path` in the `params`
.getElement <- function(params, path) {
  if (length(path) == 1) {
    params[[unlist(path)]] 
  } else {
    Recall(params[[unlist(path[1])]], path[-1])
  }
}

# alternative implementation
# get the element value located on the `path` in the `params`
.getElement2 <- function(params, path) {
  Reduce(function(z, path) z[[path]], init = params, path)
}


#' Estimate profile likelihood for gene parameters (all other fixed)
#'
#' @param pd the \link{PulseData} object
#' @param fit result output from the \link{fitModel} function
#' @param geneIndex integer(1); the row index in the count table 
#'   where the data for the given gene are located
#' @param parName a character; 
#'   the names of the gene-specific parameter (e.g. "mu")
#' @param options the options list used for the \link{fitModel} function call
#' @param ... other parameters to pass to \link{runPL},
#'  i.e. `logScale` (logical) and number of points `numPoints`
#' @inheritParams runPL
#' @return a data.frame; the column `logl` corresponds to the -log(likelihood) 
#'   function values. the other represents the profiled parameter values.
#' @export
#' @rdname profile
#'
profileGene <- function(parName,
                            geneIndex,
                            pd, 
                            par,
                            interval, 
                            options,
                            ...
                            ) {
  pL <- plGene(parName, geneIndex, par, pd, options)
  parValue <- par[[parName]][geneIndex]
  if (missing(interval))
    interval <- c(.1,2) * parValue
  result <- runPL(pL, interval, ...)
  names(result)[1] <- parName
  result
}

#' Profile
#'
#' @inheritParams profileOnlyGene
#' @inheritParams pl
#' @param parName 
#'
#' @return a data.frame; the column `logl` corresponds to the -log(likelihood) 
#'   function values. the other represents the profiled parameter values.
#' @export
#' @rdname profile
#'
profile <- function(paramPath,
                    pd,
                    par,
                    options,
                    interval,
                    namesToOptimise = names(options$lb), 
                    ...) {
  pL <- pl(paramPath, par, pd, options, namesToOptimise)
  parValue <- .getElement2(par, paramPath)
  if (missing(interval))
    interval <- c(.1,2) * parValue
  result <- runPL(pL, interval, ...)
  names(result)[1] <- as.character(paramPath[[1]]) 
  result
}

#' Compute profile likelihood on the interval
#'
#' @param pL a function (x) --> double(1) returning the -log(likelihood).
#'   For example, it can be a function returned by 
#'   \link{`pl`} or \link{`plGene`}.
#' @param interval double(2) vector; left and right boundaries to calculate the
#'   profile likelihood for the parameter given in `parName`
#' @param logScale a logical (default: `FALSE`). Should points on the 
#'   `interval` be positioned at log scale?
#' @param numPoints the number of points to position at the `interval` for
#'   profile likelihood calculations
#'
#' @return a data.frame; the first column consists of the parameter values,
#'   the second ("logL") is for -log(likelihood) values.
#' 
#' @export
#'
runPL <- function(pL,  interval, logScale = FALSE, numPoints = 21) {
  if (logScale) {
    profileParam <- exp(data.frame(x = seq(
      log(interval[1]), log(interval[2]), length.out = numPoints
    )))
  } else {
    profileParam <-
      data.frame(x = seq(interval[1], interval[2], length.out = numPoints))
  }
  res <- vapply(profileParam[,1], pL, double(1))
  profileParam$logL <- res
  profileParam
}

#' Return a profile likelihood function for further use.
#' 
#' A profile likelihood function returned by the `pl` considers parameters of
#' other genes. In contrast, the one from the `plGene` assumes
#'  parameters of other genes,  shared parameters and 
#'  the normalisation factors to be fixed. \cr
#' The profile likelihood is estimated by numerical optimisation, and it can
#' be sensitive to the initial values. The optimisation is repeated 
#' `options$replicate` times (default: 10), by adding a random values to the 
#' the optimum  parameters at a scale specified in the options$jitter.
#' 
#' 
#' @param paramPath a list with names and indexes in order to locate the 
#' parameter in the `par` argument (e.g. `list("mu", 1)` corresponds to
#' the "mu" parameter value for the first gene, i.e.  
#' `\verb{par[["mu"]][[1]]}` 
#' @param parName a character, e.g. "mu"
#' @param geneIndex an integer(1); a row index which corresponds to 
#'   the investigating gene
#' @param par a result of the \link{fitModel} function
#' @param pd a \link{PulseData} object
#' @param options an option list used for the \link{fitModel} call;
#'   additional options can be specified:
#'  - jitter, a double (default: .1)
#'  - replicate, an integer number of repeating the optimisation from random 
#'    points
#'  - absolute, a logical (default: FALSE); if FALSE, the likelihood value at
#'  the optimal point (i.e. `par`) is substracted from the returned value.
#'  In this case, the value of the returned function at the point `par` is 0.
#' @param freeParams which parameters are optimised (i.e. not fixed);
#'   by default they are derived from the names of the boundaries
#'
#' @details The randomisation of the parameter x is made as 
#' $ x^* = (1 + a) * x$, where $a$ is a random value from the uniform
#' distribution (0,`options$jitter`). 
#'
#' @return a function with the calling convention as   
#'   `f(x = double(1)) --> double(1))`,  
#'   which return the value of -log(likelihood(mu = x))
#' 
#' @name profiles
NULL

.addDefaultsPL <- function(options) {
  defaults <- list(
    replicates = 10,
    jitter = .1,
    absolute = FALSE)
  notSet <- setdiff(names(defaults),names(options))
  options[notSet] <- defaults[notSet]
  options
}

#' @export
#' @rdname profiles
plGene <- function(parName,
                   geneIndex,
                   par,
                   pd,
                   options) {
  options <- .addDefaultsPL(options)
  knownNames <- .getKnownNames(par, options)
  namesToOptimise <- setdiff(.getGeneToFitNames(par, knownNames), parName)
  options <- normaliseBoundaries(options, par, pd)
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(options$lb[namesToOptimise])
  ub <- as.data.frame(options$ub[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  knownNames <- c(knownNames, parName)
  fixedPars <- par[!names(par) %in% namesToOptimise]
  fixedPars[knownNames] <- lapply(par[knownNames], `[[`, geneIndex)
  initValues <- unlist(lapply(par[namesToOptimise], `[[`, geneIndex))
  optimum <- ifelse(options$absolute,
                    0,
                    objective(initValues, pd$counts[geneIndex, ], fixedPars))
  function(x) {
    fixedPars[parName] <- x
    min(replicate(options$replicates,{
      jitterCoeffs <- 1 + runif(length(initValues), 
                                -options$jitter, options$jitter)
    .fitGene(initValues * jitterCoeffs, geneIndex, objective, lb, ub, 
             fixedPars, pd$counts)$value - optimum
    })) }
}


#' @export
#' @rdname profiles
pl <- function(paramPath,
               par,
               pd,
               options,
               freeParams = names(options$lb)) {
  options <- .addDefaultsPL(options)
  knownNames <- .getKnownNames(par, options)
  options <- normaliseBoundaries(options, par, pd)
  optimalValue <- .getElement2(par, paramPath)
  boundaries <- lapply(c(lb = "lb", ub = "ub"), function(side) {
    b <- options[[side]][freeParams]
    unlist(b)[.getFreeIndexes(b, paramPath, freeParams)]
  })
  freeInd <- .getFreeIndexes(par, paramPath, freeParams)
  objective <- function(x, params) {
    p <- unlist(params[freeParams])
    p[freeInd] <- x
    params[freeParams] <- relist(p, params[freeParams])
    -evaluateLikelihood(params, pd)
  }
  optimisationStart <- unlist(par[freeParams])[freeInd]
  optimum <- -evaluateLikelihood(par, pd)
  function(x) {
    par <- .assignElement(par, paramPath, x)
    stats::optim(
      optimisationStart,
      objective,
      method = "L-BFGS-B",
      control = list(parscale = optimisationStart),
      lower = boundaries$lb,
      upper = boundaries$ub,
      params = par
    )$value - optimum
  }
}

# fit params for i-th gene
# initPars is a vector of double
# objective is a function to optimise 
# the calling convention if f(x, counts, fixedPars),
# where x are the parameters to fit, fixed is a character vector of gene-
# sepecific parameters which are fixed
.fitGene <- function(initPars, i, objective, lb, ub, fixedPars, counts) {
  stats::optim(
    initPars,
    objective,
    method = "L-BFGS-B",
    control = list(parscale =  initPars), 
    lower = lb[i, ],
    upper = ub[i, ],
    counts = counts[i, ],
    fixedPars = fixedPars
  )
}

# returns indexes of the parameters in the unlist(par[namesToOptimise])
# additionally, the parameter at paramPath is excluded if it is present in
# par[namesToOptimise]
.getFreeIndexes <- function(par, paramPath, namesToOptimise) {
  freePars <- par[namesToOptimise]
  fixedIndexes <- .getIndex(paramPath, freePars)
  freeInd <- unlist(relist(seq_along(unlist(freePars)),freePars))
  freeInd <- freeInd[freeInd != fixedIndexes]
  if (!is.null(par$normFactors)) {
    freeInd <- freeInd[freeInd != .getIndex(list("normFactors", 1, 1), freePars)]
  }
  freeInd
}


#' Plot the profile likeliihood
#'
#' @param pl  a result from the \link{profileOnlyGene} frunction
#' @param confidence a confidence level for the likelihood threshold line
#'   (default .95)
#'
#' @return used for its side effect
#' @export
#' @rdname profile
#'
plotPL <- function(pl, confidence=.95, ...) {
  parName <- names(pl)[which(names(pl) != "logL")]
  par(mar = c(4,4,1,1))
  plot(
    x = pl[[parName]],
    y = pl[["logL"]],
    type = "l",
    xlab = "", #parName,
    ylab =  "", #"-logL",
    bty = "l",
    ...
  )
  mtext(parName,1,2, family = "serif")
  mtext("-logL",2,2.5, family = "serif")
  abline(h = qchisq(confidence,1)/2, col= rgb(1,.2,0), lwd=2)
}


#' Estimate confidence intervals
#'
#' @inheritParams plGene
#' @inheritParams getCI
#' @param geneIndexes a vector;  corresponds to the indexes of genes, for which the
#' confidence intervals must be computed
#' @details
#' @return 
#' - `ciGene`: 
#'   a data.frame with two columns (left, right) confidence boundaries
#'   in the order of geneIndexes;  
#' - `ci`: a vector of two numbers;
#' @export
#' @rdname confint
#'
ciGene <- function(parName, geneIndexes, pd,  par, options, interval,
                   confidence = .95) {
  options <- normaliseBoundaries(options, par, pd)
  if (missing(interval)) {
    interval <- cbind(options$lb[[parName]], options$ub[[parName]])
  }
  interval_ <- interval
  if (!is.null(options$parallel) && options$parallel) {
    applyfun <- mclapply
  } else {
    applyfun <- lapply
  }
  result <- lapply(geneIndexes, function(geneIndex) {
    if (!is.null(dim(interval)[1]))
      interval_ <- interval[geneIndex,]
    pL <- plGene(parName, geneIndex, par, pd, options)
    getCI(
      pL = pL,
      optimum = par[[parName]][geneIndex],
      interval = interval_,
      confidence = confidence
    )
  })
  do.call(rbind, result)
}

#' @export
#' @rdname confint
ci <- function(paramPath, pd, par, options, freeParams,
               interval, confidence = .95) {
  options <- normaliseBoundaries(options, par, pd)
  if (missing(interval)) {
    interval <- cbind(.getElement2(options$lb, paramPath),
                      .getElement2(options$ub, paramPath))
  }
  if (missing(freeParams)) {
    pL <- pl(paramPath, par, pd, options)
  } else {
    pL <- pl(paramPath, par, pd, options, freeParams)
  }
  getCI(
    pL = pL,
    optimum = .getElement2(par, paramPath),
    interval = interval,
    confidence = confidence
  )
}

#'  Estimate confidence interval
#'
#' @param pL a profile likelihood function
#' @param optimum an optimal value of the parameter (e.g. from fitting)
#' @param confidence confidence level for interval estimation
#' @param interval a vector of too numbers, define (min, max) allowed 
#'   parameter values
#'
#' @return a vector of two: c(min, max);
#' @details  if the likelihood does not reach 
#'   the needed level (i.e. the CI exceeds the limitations in the `interval`,
#'   NA is returned. Hence, in case of both ends, the return value is c(NA, NA).
#' @keywords internal
#' 
getCI <- function(pL, optimum, confidence, interval) {
  threshold <- qchisq(confidence, 1)/2
  objective <- function(x) pL(x) - threshold
  optimalObjective <- objective(optimum)
  ci <- c(NA, NA)
  if (optimalObjective * objective(interval[1]) < 0)
    ci[1] <- uniroot(objective, c(interval[1], optimum))$root
  if (optimalObjective * objective(interval[2]) < 0)
    ci[2] <- uniroot(objective, c(optimum, interval[2]))$root
  ci
}


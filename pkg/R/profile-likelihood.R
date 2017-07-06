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
#'
profileOnlyGene <- function(pd, par,
                            geneIndex,
                            parName,
                            options,
                            interval, 
                            ...
                            ) {
  pL <- plFixed(parName, par, options, pd, geneIndex)
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
#'
profile <- function(paramPath,
                    pd,
                    par,
                    options,
                    interval,
                    parName = as.character(paramPath[[1]]) ,
                    namesToOptimise = names(options$lb), 
                    ...) {
  pL <- pl(paramPath, par, options, pd,namesToOptimise )
  result <- runPL(pL, interval, ...)
  names(result)[1] <- parName
  result
}

#' Compute profile likelihood on the interval
#'
#' @param pL a function (x) --> double(1) returning the -log(likelihood)
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

#' Get the profile likelihood function (all other parameter fixed)
#'
#' @param parName a character, e.g. "mu"
#' @param par a result of the \link{fitModel} function
#' @param options an option list used for the \link{fitModel} call
#' @param pd a \link{PulseData} object
#' @param geneIndex an integer(1); a row index which corresponds to 
#'   the investigating gene
#'
#' @return a function with the calling convention as   
#'   `f(x = double(1)) --> double(1)`,  
#'   which return the value of -log(likelihood(mu = x))
#' @export
#'
plFixed <- function(parName,
                    par,
                    options,
                    pd,
                    geneIndex) {
  knownNames <- .getKnownNames(par, options)
  namesToOptimise <- setdiff(.getGeneToFitNames(par, knownNames), parName)
  .pLfunction(options, par, pd, parName, geneIndex, namesToOptimise, knownNames)
}

#' Get the profile likelihood function (consider the parameters of other genes
#' as not fixed)
#'
#' @param paramPath a list with names and indexes in order to locate the 
#' parameter of the profile
#' @param par a result of the \link{fitModel} function
#' @param options an option list used for the \link{fitModel} call
#' @param pd a \link{PulseData} object
#' @param namesToOptimise which parameters are optimised (i.e. not fixed);
#'   by default they are derived from the names of the boundaries
#'
#' @return a function with the calling convention as   
#'   `f(x = double(1)) --> double(1)`,  
#'   which return the value of -log(likelihood(mu = x))
#'
#' @details `paramPath` is defined as the sequence of selections performed on 
#' the parameter list (`par`). For example, list("mu", 1) will point to the 
#' "mu" parameter of the first gene. 
#' 
#' @export
#'
pl <- function(paramPath,
               par,
               options,
               pd,
               namesToOptimise = names(options$lb)) {
  knownNames <- .getKnownNames(par, options)
  #.pLfunctionTotal(options, par, pd, paramPath, namesToOptimise, knownNames)
  options <- normaliseBoundaries(options, par, pd)
  optimalValue <- .getElement2(par, paramPath)
  boundaries <- lapply(c("lb", "ub"), function(side) {
    b <- options[[side]][namesToOptimise]
    unlist(b)[.getFreeIndexes(b, paramPath, namesToOptimise)]
  })
  names(boundaries) <- c("lb", "ub")
  freeInd <- .getFreeIndexes(par, paramPath, namesToOptimise)
  objective <- function(freeParams, params) {
    x <- unlist(params[namesToOptimise])
    x[freeInd] <- freeParams
    params[namesToOptimise] <- relist(x, params[namesToOptimise])
    -evaluateLikelihood(params, pd)
  }
  optimisationStart <- unlist(par[namesToOptimise])[freeInd]
  optimum <- -evaluateLikelihood(par, pd)
  function(x) {
    par <- .assignElement(par, paramPath, x)
    stats::optim(
      optimisationStart,
      objective,
      method = "L-BFGS-B",
      control = list(parscale = (boundaries$lb + boundaries$ub) / 2),
      lower = boundaries$lb,
      upper = boundaries$ub,
      params = par
    )$value - optimum
  }
  
}

# fit params for i-th gene
# p is a data.frame with the being fitted parameters by column
# objective is a function to optimise 
# the calling convention if f(x, counts, fixedPars),
# where x are the parameters to fit, fixed is a character vector of gene-
# sepecific parameters which are fixed
.fitGene <- function(p, i, objective, lb, ub, fixedPars, counts) {
  stats::optim(
    unlist(p[i,]),
    objective,
    method = "L-BFGS-B",
    control = list(parscale = p[i, ]),
    lower = lb[i, ],
    upper = ub[i, ],
    counts = counts[i, ],
    fixedPars = fixedPars
  )
}

.pLfunction <- function(options,
                       par,
                       pd,
                       paramName,
                       geneIndex,
                       namesToOptimise,
                       knownNames) {
  options <- normaliseBoundaries(options, par, pd)
  optimalValue <- par[[paramName]][geneIndex]
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(options$lb[namesToOptimise])
  ub <- as.data.frame(options$ub[namesToOptimise])
  initValues <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  par[namesToOptimise] <- NULL
  knownNames <- c(knownNames,paramName)
  fixedPars <- par
  fixedPars[knownNames] <- lapply(par[knownNames], `[[`, geneIndex)
  fixedPars[paramName] <- optimalValue
  optimum <- .fitGene(initValues, geneIndex, objective, lb, ub, 
             fixedPars, pd$counts)$value
  function(x) {
    fixedPars[paramName] <- x
    .fitGene(initValues, geneIndex, objective, lb, ub, 
             fixedPars, pd$counts)$value - optimum
  }
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


#' Estimate CI when parameters of other genes are fixed
#'
#' @inheritParams plFixed
#' @inheritParams getCI
#' @param geneIndexes an integer vector; for which genes subset confidence
#'   intervals are to be calculated
#' @param interval 
#'
#' @return a data.frame with two columns (left, right) confidence boundaries
#'   in the order of geneIndexes
#' @export
#'
estimateCIFixed <- function(pd,
                            fit,
                            geneIndexes,
                            parName,
                            options,
                            confidence=.95,
                            interval) {
  if (missing(interval)) {
    interval <- cbind(options$lb[[parName]], options$ub[[parName]])
  }
  result <- lapply(geneIndexes, function(geneIndex) {
    pL <- plFixed(
      parName   = parName,
      par       = fit,
      options   = options,
      pd        = pd,
      geneIndex = geneIndex
    )
    if (!is.null(dim(interval)) && (dim(interval) > 1)) {
      interval <- interval[geneIndex, ]
    }
    getCI(
      pL = pL,
      optimum = fit[[parName]][geneIndex],
      interval = interval,
      confidence = confidence
    )
  })
  do.call(rbind, result)
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


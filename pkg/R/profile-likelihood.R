.getIndex <- function(path, example) {
  x <- unlist(example)
  x[] <- seq_along(x)
  index <- .getElement(relist(x, example), path)
  if (length(index) > 1) 
    stop("Not complete path to the element (length of selection > 1)")
  index
}
  
.assignElement <- function(x, i, value) {
  if (length(i) == 1) {
    x[[unlist(i)]] <- value 
    x
  } else {
    x[[unlist(i[1])]] <- Recall(x[[unlist(i[1])]], i[-1], value)
    x 
  }
}

.getElement <- function(x, i) {
  if (length(i) == 1) {
    x[[unlist(i)]] 
  } else {
    Recall(x[[unlist(i[1])]], i[-1])
  }
}

.getElement2 <- function(x,i) {
  Reduce(function(z, i) z[[i]], init = x, i)
}


#' Profile likelihood for a single parameter
#'
#' @param pd a PulseData object
#' @param fit a fitting result
#' @param opts options object 
#' @param var a list, describing position of the parameter in `fit`
#' @param interval a vector of two numbers. Defines boundaries for profile
#' likelihood estimation
#'
#' @return a data.frame ["value", "logL"]
#'
profile <- function(pd, result, opts, var, interval, N = 20) {
  x <- seq(interval[1], interval[2], seq.length = N) 
  logL <- lapply(x, function(value) {
    .profileOptim(path = var,
                  value = value,
                  opts = opts)
  }) 
  data.frame(values = values, logL = logL)
}

#' Estimate profile likelihood for gene parameters (all other fixed)
#'
#' @param pd the \link{`PulseData`} object
#' @param fit result output from the \link{`fitModel`} function
#' @param geneIndex integer(1); the row index in the count table 
#'   where the data for the given gene are located
#' @param parName a character; 
#'   the names of the gene-specific parameter (e.g. "mu")
#' @param options the options list used for the \link{`fitModel`} function call
#' @param interval double(2) vector; left and right boundaries to calculate the
#'   profile likelihood for the parameter given in `parName`
#' @param logScale a logical (default: `FALSE`). Should points on the 
#'   `interval` be positioned at log scale?
#' @param numPoints the number of points to position at the `interval` for
#'   profile likelihood calculations
#'
#' @return a data.frame; the column `logL` corresponds to the -log(likelihood) 
#'   function values. The other represents the profiled parameter values.
#' @export
#'
profileOnlyGene <- function(pd,
                            par,
                            geneIndex,
                            parName,
                            options,
                            interval,
                            logScale = FALSE,
                            numPoints = 21) {
  if (logScale) {
    profileParam <- exp(data.frame(x = seq(
      log(interval[1]), log(interval[2]), length.out = numPoints
    )))
  } else {
    profileParam <-
      data.frame(x = seq(interval[1], interval[2], length.out = numPoints))
  }
  names(profileParam) <- parName
  pL <- plFixed(parName, par, options, pd, geneIndex)
  res <- vapply(profileParam[,1], pL, double(1))
  profileParam$logL <- res
  profileParam
}

#' Get the profile likelihood function (all other parameter fixed)
#'
#' @param parName a character, e.g. "mu"
#' @param par a result of the \link{`fitModel`} function
#' @param options an option list used for the \link{`fitModel`} call
#' @param pd a \link{`PulseData`} object
#' @param geneIndex an integer(1); a row index which corresponds to 
#'   the investigating gene
#'
#' @return a function with the calling convention as   
#'   `f(x = double(1)) --> double(1)`,  
#'   which return the value of -log(likelihood(mu = x))
#' @export
#'
#' @examples
plFixed <- function(parName,
                    par,
                    options,
                    pd,
                    geneIndex) {
  knownNames <- .getKnownNames(par, options)
  namesToOptimise <- setdiff(.getGeneToFitNames(par, knownNames), parName)
  pLfunction(options, par, pd, parName, geneIndex, namesToOptimise, knownNames)
}

pLfunction <- function(options,
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
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  par[namesToOptimise] <- NULL
  fixedPars <- par
  knownNames <- c(knownNames,paramName)
  fixedPars[knownNames] <-
    lapply(par[knownNames], `[[`, geneIndex)
  fixedPars[paramName] <- optimalValue
  optimum <-
    .fitGene(p, geneIndex, objective, lb, ub, fixedPars, pd$counts)$value
  function(x) {
    fixedPars[paramName] <- x
    .fitGene(p, geneIndex, objective, lb, ub, fixedPars, pd$counts)$value - optimum
  }
}

#' Plot the profile likeliihood
#'
#' @param pl  a result from the \link{`profileOnlyGene`} frunction
#' @param confidence a confidence level for the likelihood threshold line
#'   (default .95)
#'
#' @return used for its side effect
#' @export
#'
plotPL <- function(pl, confidence=.95) {
  parName <- names(pl)[which(names(pl) != "logL")]
  par(mar = c(4,4,1,1))
  plot(
    x = pl[[parName]],
    y = pl[["logL"]],
    type = "l",
    xlab = "", #parName,
    ylab =  "", #"-logL",
    bty = "l"
  )
  mtext(parName,1,2, family = "serif")
  mtext("-logL",2,2.5, family = "serif")
  abline(h = qchisq(confidence,1)/2, col= rgb(1,.2,0), lwd=2)
}


#' Estimate CI when parameters of other genes are fixed
#'
#' @param pd 
#' @param fit 
#' @param geneIndexes 
#' @param parName 
#' @param options 
#' @param confidence 
#' @param interval 
#'
#' @return
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
    .getCI(
      pL = pL,
      optimum = fit[[parName]][geneIndex],
      interval = interval,
      confidence = confidence
    )
  })
  do.call(rbind, result)
}

.getCI <- function(pL, optimum, confidence, interval) {
  threshold <- qchisq(confidence, 1)/2
  ciLeft <- uniroot(function(x) pL(x) - threshold, c(interval[1], optimum))$root
  ciRight <- uniroot(function(x) pL(x) - threshold, c(optimum, interval[2]))$root
  c(ciLeft, ciRight)
}


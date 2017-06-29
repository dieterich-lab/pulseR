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
#' @export
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

if (fixNorms) {
  p <- fitParamsSeparately(
    pd = pd,
    par = par,
    knownNames = known,
    namesToOptimise = toFit,
    options = opts,
    indexes = geneIndex
  )
  
} else {
  
}

profileOnlyGene <- function(pd,
                            fit,
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
  knownNames <- .getKnownNames(fit, options)
  namesToOptimise <- setdiff(.getGeneToFitNames(fit, knownNames), parName)
  .profileGene(pd,
               fit,
               knownNames,
               namesToOptimise,
               profileParam,
               options,
               geneIndex) 
}

.profileGene <- function(pd,
                        par,
                        knownNames,
                        namesToOptimise,
                        profileParam,
                        options,
                        geneIndex) {
  options <- normaliseBoundaries(options, par, pd)
  optimalValue <- par[[names(profileParam)]][geneIndex]
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(options$lb[namesToOptimise])
  ub <- as.data.frame(options$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  par[namesToOptimise] <- NULL
  fixedPars <- par
  knownNames <- c(knownNames, names(profileParam))
  fixedPars[knownNames] <- lapply(par[knownNames], `[[`, geneIndex)
  fixedPars[names(profileParam)] <- optimalValue
  optimum <-
    .fitGene(p, geneIndex, objective, lb, ub, fixedPars, pd$counts)$value
  pL <- double(length(profileParam[, 1]))
  for (i in seq_along(pL)) {
    fixedPars[names(profileParam)] <- profileParam[i, , drop = FALSE]
    pL[i] <-
      .fitGene(p, geneIndex, objective, lb, ub, fixedPars, pd$counts)$value
  }
  pL <- pL - optimum
  profileParam$logL <- pL
  profileParam
}

fit <- result
geneIndex <- 10
parName <- "b"
interval <- c(0.9,1.1) * fit[[parName]][geneIndex]

pl <- profileOnlyGene(
  pd,
  fit,
  geneIndex, parName, options, interval, numPoints = 201)
plot(pl, type="l")
abline(h = qchisq(.95,1)/2, col=5)

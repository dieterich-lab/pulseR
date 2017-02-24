
#' Create an object for pulse-change count data
#'
#' @param counts a matrix; column names correspond to sample names
#' @param conditions a data.frame;
#'   the first column corresponds to the conditions given in \code{formulas}.
#' @param formulas a list, created by \code{\link{MeanFormulas}}
#' @param formulaIndexes list of lists; defines indexes of formulas 
#' used for calculation of the expected read number
#' @param spikeins 
#' @param fractions a formula, e.g. ~ condition + time (if spike-ins are
#' not provided).
#' the names used in the \code{fractions} defines different fractions,
#' which should have distinct coefficients for mean expression fitting.
#'
#' @return an object of class "PulseData"
#' @export
#'
PulseData <- function(counts,
                      conditions,
                      formulas,
                      formulaIndexes = NULL,
                      spikeins = NULL,
                      groups = NULL) {
  e <- new.env()
  if (is.null(formulaIndexes))
    formulaIndexes <- match(conditions[, 1], names(formulas))
  class(e) <- "PulseData"
  e$user_conditions <- conditions
  e$counts <- as.matrix(counts)
  known <- addKnownToFormulas(formulas, formulaIndexes, conditions)
  e$conditions <- conditions
  e$formulas <- known$formulas
  e$formulaIndexes <- known$formulaIndexes
  e$user_formulas <- formulas
  if (!is.null(groups)) {
    if (is(groups, "formula"))
      groups <- interaction(conditions[, all.vars(groups)])
    g <- makeGroups(e, groups)
    e$normCoeffs <- g$normCoeffs
    e$normCoeffIndexes <- g$normCoeffIndexes
  }
  #normalise(e)
  e
}



#' @export
print.PulseData <- function(x,...){
  cat("PulseData object")
}

#' Calculate normalisation factors for columns in a matrix
#'
#' @param count_data a matrix; columns correspond to samples.
#'
#' @return a vector of doubles with normalisation factors ordered as columns in 
#'   the \code{count_data}
#'
#' @importFrom  stats median
findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    finitePositive <- is.finite(loggeomeans) & x > 0
    if (any(finitePositive))
      res <- exp(median((log(x) - loggeomeans)[finitePositive], na.rm = TRUE))
    else {
      print(count_data[1:6,])
      stop("Can't normalise accross a condition. 
              Too many zero expressed genes. ")
    }
    res  
  })
  deseqFactors
}

#' Calculate normalisation factors
#'
#' @param count_data integer matrix, colnames correspond to samples 
#'   (rownames in \code{conditions})
#' @param conditions factors to split samples for normalisation
#' @return vector of double; normalisation factors in the same order as 
#'   columns in the \code{count_data}
findDeseqFactorsForFractions <- function(count_data, conditions) {
    deseqFactors <- lapply(
      split(colnames(count_data), conditions),
      function(samples) {
        findDeseqFactorsSingle(count_data[, samples, drop = FALSE])
      })
    names(deseqFactors) <- NULL
    unlist(deseqFactors)[colnames(count_data)]
}

# evaluate formulas in the environment of known params from the conditions
# Returns list( [evaluated formulas], [conditions as vector])

addKnownToFormulas <- function(forms, formulaIndexes, conditions) {
  uc <- unique(conditions)
  newIndexes <- list()
  newForms <- list()
  for (i in seq_along(uc[, 1])) {
    f <- forms[formulaIndexes[[as.character(uc[i, 1])]]]
    newNames <- as.character(
      interaction(c(list(names(f)), uc[i, -1])))
    res <- lapply(f, substitute_q, env = as.list(uc[i, ]))
    names(res) <- newNames
    newForms[names(res)] <- res
    newIndexes[[i]] <- newNames
  }
  names(newIndexes) <- interaction(uc)
  newIndexes <- multiplyList(newIndexes, interaction(conditions))
  newIndexes <- names2numbers(newIndexes, names(newForms))
  list(formulas = newForms, formulaIndexes = newIndexes)
}

makeGroups <- function(pd, normGroups) {
  normCoeffs <- pd$formulaIndexes[match(unique(normGroups), normGroups)]
  names(normCoeffs) <- unique(normGroups)
  normCoeffs <-
    relist(seq_along(unlist(normCoeffs)), normCoeffs)
  normCoeffIndexes <- multiplyList(normCoeffs, normGroups)
  list(normCoeffs = normCoeffs,
       normCoeffIndexes = normCoeffIndexes)
}

names2numbers <- function(nameLists, nameVector){
  lapply(nameLists, match, nameVector)
}

#' Create a test count data
#'
#' @param formulas a list
#' @param formulaIndexes list of lists; defines indexes of formulas 
#' used for calculation of the expected read number
#' @param normFactors list of vectors; normalisation factors, if
#' known
#' @param par a list of named parameters; gene-specific parameters 
#' are vectors
#' @param conditions a condition data.frame
#' @return matrix of counts with the order of columns as in conditions 
#' @importFrom stats rnbinom
#' @export
#'
generateTestDataFrom <- function(formulas,
                                 formulaIndexes,
                                 normFactors,
                                 par,
                                 conditions) {
  if (all(vapply(formulaIndexes, is.character, logical(1))))
    formulaIndexes <- names2numbers(formulaIndexes, names(formulas))
  known <- addKnownToFormulas(formulas, formulaIndexes, conditions)
  formulas <- known$formulas
  indexes <- known$formulaIndexes
  evaled <- lapply(formulas, eval, env=par) 
  means <- sample_means(evaled, indexes, normFactors)
  counts <- matrix(rnbinom(length(means), mu = means, size = par$size),
         ncol = length(conditions[,1]))
  counts
}

multiplyList <- function(source, pattern) {
  res <- list()
  for (i in seq_along(pattern)) {
    res[[i]] <- source[[as.character(pattern[i])]]
  }
  names(res) <- pattern
  res
}

shrinkList <- function(l){
 l[unique(names(l))] 
}
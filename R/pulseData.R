
#' Create an object for pulse-change count data
#'
#' @param counts a matrix; column names correspond to sample names
#' @param conditions a data.frame;
#'   the first column corresponds to the conditions given in \code{formulas}.
#' @param formulas a list, created by \code{\link{MeanFormulas}}
#' @param formulaIndexes list of lists; defines indexes of formulas 
#' used for calculation of the expected read number
#' @param spikeins 
#' @param fractions a formula or a vector, e.g. ~ condition + time (if spike-ins are
#' not provided). The vector length must be the same as the sample number.
#' the names used in the \code{fractions} defines different fractions,
#' which should have distinct coefficients for mean expression fitting.
#'
#' @return an object of class "PulseData"
#' @export
#' 
#' @examples 
#' 
#' formulaIndexes <- list(
#'   total_fraction = 'total',
#'   flow_through   = c('unlabelled', 'labelled'),
#'   pull_down      = c('labelled', 'unlabelled'))
#'   
#' # spike-ins set up:
#' 
#' refGroup <- "total_fraction"
#' 
#' labelled <- c("spike1", "spike2") 
#' unlabelled <- c("spike3", "spike4") 
#' 
#' spikeLists <- list(
#' #  total samples are normalised using all spike-ins
#'   total_fraction = list(c(unlabelled, labelled)),
#' # for every item in formulaIndexes we have a set of spike-ins:   
#'   flow_through   = list(unlabelled, labelled),
#'   pull_down      = list(labelled, unlabelled))
#'   
#' # argument for the function: 
#' spikeins <- list(refGroup = refGroup,
#'                  spikeLists = spikeLists)
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
  # generate partially evaluated formulas for every condition
  # and create new formula indexes for every of this conditions
  known <- addKnownToFormulas(formulas, formulaIndexes, conditions)
  e$conditions <- conditions
  e$rawFormulas <- known$formulas
  # the compiled formulas are used during likelihood calculation
  e$formulas <- lapply(known$formulas,
                     compiler::compile,
                     options = list(suppressAll = TRUE))
  e$formulaIndexes <- known$formulaIndexes
  e$user_formulas <- formulas
  # create a skeleton for the sequencing depth normalisation coefficients
  e$depthNormalisation <- assignList(e$formulaIndexes, 1)
  if (!is.null(spikeins)) {
    # normalise to spike-ins and remove them from the count table
    e$depthNormalisation <- normaliseWithSpikeIns(
      e, spikeins$refGroup, spikeins$spikeLists)
    refSpikes <- unlist(spikeins$spikeLists[[spikeins$refGroup]])
    if (is.character(refSpikes))
      refSpikes <- which(rownames(e$counts) %in% refSpikes)
    e$counts <- e$counts[-refSpikes,]
  } else {
    if (is.null(groups))
      groups <- seq_along(conditions[,1])
    if (is(groups, "formula"))
      groups <- interaction(conditions[, all.vars(groups)])
    e$groups <- groups
    g <- makeGroups(e, groups)
    e$interSampleCoeffs  <- g$normCoeffs
    e$interSampleIndexes <- g$normCoeffIndexes
    deseqFactors <- normaliseNoSpikeins(e, groups)
    for (i in seq_along(e$depthNormalisation)) {
      e$depthNormalisation[[i]][] <- deseqFactors[i]
    }
  }
  e
}


#' @export
print.PulseData <- function(x,...){
  cat("PulseData object \n")
  cat(paste0(dim(x$counts)[1], " genes X ", 
             dim(x$counts)[2], " samples "))
 print(table(conditions[,1]))
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
  deseqFactors <- apply(count_data, 2, deseq, loggeomeans = loggeomeans) 
  deseqFactors
}

deseq <- function(x, loggeomeans) {
  finitePositive <- is.finite(loggeomeans) & x > 0
  if (any(finitePositive)) {
    res <- exp(median((log(x) - loggeomeans)[finitePositive], na.rm = TRUE))
  } else {
    print(head(x))
    stop("Can't normalise accross a condition.
         Too many zero expressed genes. ")
  }
  res
}

normaliseWithSpikeIns <- function(pd, refGroup, spikeLists){
  refSpikes <- unlist(spikeLists[[refGroup]])
  spikeCounts <- pd$counts[refSpikes,]
  spikeLists <- lapply(spikeLists,
                       function(l) {
                         lapply(l, function(spikes)
                           spikes == refSpikes)
                       })
                       
  superSample <- rowMeans(
    log(pd$counts[refSpikes, pd$conditions$condition == refGroup]))
  lapply(seq_along(pd$conditions[, 1]),
         function(i) {
           sampleSpikes <- spikeLists[[as.character(pd$conditions[i, 1])]]
           vapply(sampleSpikes, function(spikes) {
             deseq(spikeCounts[spikes, i], superSample)
           }, double(1))
         })
}


#' Performs sequencing depth normalisation using the DESeq procedure.
#'
#' @param pd a \code{\link{PulseData}} object 
#' @param groups a vector for splitting objects to groups
#'
#' @return a vector with the coefficient for every sample
#'
normaliseNoSpikeins <- function(pd, groups){
  factors <- double(length(groups))
  for (g in unique(groups)) {
    factors[groups == g] <- findDeseqFactorsSingle(
      pd$counts[, groups == g, drop = FALSE])
  }
  factors
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

#' Evaluate formulas in the environment of known params from the conditions
#' 
#' If, for example, a labelled fraction is estimated at several time points
#' (0hr, 2hr, 4hr), corresponding partially evaluated formulas will be
#' created. In this case, a time variable (e.g. "t") will be substituted by
#' its values from the corresponding column in the conditions data.frame.
#' Hence several partially evaluated formulas will be created for every
#' combination of variable values described in the conditions data.frame.
#'
#' @param formulas a named list with unevaluated expressions
#'  for expression levels 
#'  (e.g. describing amounts of "labelled", "total" etc. RNA)
#' @param formulaIndexes a list describing which formulas to use for
#' mean read number calculation for fractions defined by the names of list items
#' @param conditions a data.frame with the first column corresponding to 
#' names in formulasIndexes (e.g. "total", "pull_down")
#'
#' @return a list with two items:
#'  - list of partially evaluated formulas
#'  - a vector of conditions generated from combination of columns 
#' @export
#'
#' @examples
#' 
#' formulas <- MeanFormulas(total = m, label = m * exp(-d*t))
#' formulaIndexes <- list(
#'   total = 'total',
#'   pull_down = 'label'
#' )
#' conditions <- data.frame(
#'   type = c('total', 'pull_down', 'pull_down'),
#'   t = c(0, 1, 5)
#' )
#' result <- addKnownToFormulas(formulas, formulaIndexes, conditions)
#' str(result)
#' # List of 2
#' # $ formulas      :List of 3
#' # ..$ total.0: symbol m
#' # ..$ label.1: language m * exp(-d * 1)
#' # ..$ label.5: language m * exp(-d * 5)
#' # $ formulaIndexes:List of 3
#' # ..$ total.0    : int 1
#' # ..$ pull_down.1: int 2
#' # ..$ pull_down.5: int 3
#'
addKnownToFormulas <- function(formulas, formulaIndexes, conditions) {
  uc <- unique(conditions)
  newIndexes <- list()
  newForms <- list()
  for (i in seq_along(uc[, 1])) {
    # take formulas for the current condition i
    f <- formulas[formulaIndexes[[as.character(uc[i, 1])]]]
    # a combined formula name is a concatenation with the i-th row values 
    newNames <- as.character(
      interaction(c(list(names(f)), uc[i, -1])))
    # partial evaluation in the environment of the i-th row values
    res <- lapply(f, substitute_q, env = as.list(uc[i, ]))
    names(res) <- newNames
    newForms[names(res)] <- res
    newIndexes[[i]] <- newNames
  }
  names(newIndexes) <- interaction(uc)
  # indexes must be provided for every sample
  newIndexes <- multiplyList(newIndexes, interaction(conditions))
  newIndexes <- names2numbers(newIndexes, names(newForms))
  list(formulas = newForms, formulaIndexes = newIndexes)
}

#' Constructs coefficients list and their indexes according to the grouping rule
#' 
#' It is assumed that all samples in a group belong to the same condition, 
#' i.e. the same equations are used to estimate expression levels.
#' However, they may correspond to different time points etc.
#' 
#' The returned normalisation coefficients are assign to 1.
#'
#' @param pd a \code{\link{PulseData}} object
#' @param normGroups a vector defining a rule for splitting the samples in `pd`
#'
#' @return a named list with items "normCoeffs" and "normCoeffIndexes".
#' normCoeffs is the list to store the values of the normalisation coefficients
#' for every group. 
#' normCoeffIndexes stores indexes of coefficients from unlist(normCoeffs) 
#' sample-wise, i.e. length(normCoeffIndexes) is the number of samples.
#'
makeGroups <- function(pd, normGroups) {
  # generate a list of normalisation coefficients with a proper structure
  normCoeffs <- pd$formulaIndexes[match(unique(normGroups), normGroups)]
  names(normCoeffs) <- unique(normGroups)
  # all the normalisation coefficients are numbered according 
  # to their appearance in the flatten list `unlist(normCoeffs)`
  normCoeffs <- relist(seq_along(unlist(normCoeffs)), normCoeffs)
  # normCoeffIndexes contains items for every sample
  normCoeffIndexes <- multiplyList(normCoeffs, normGroups)
  normCoeffs <- assignList(normCoeffs, 1)
  list(normCoeffs       = normCoeffs,
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
  evaled <- do.call(cbind, lapply(formulas, eval, env=par))
  pd <- list(formulas = formulas, formulaIndexes = indexes,
             depthNormalisation = normFactors)
  norms <- getNorms(pd)
  means <- sample_means(evaled, norms)
  counts <- matrix(rnbinom(length(means), mu = means, size = par$size),
         ncol = length(conditions[,1]))
  counts
}

#' A helper to generate named lists
#'
#' @param source a named lists with the values to multiply
#' @param pattern a vector of names defining how to construct a new
#' list (typically, a longer one with replicates).
#'
#' @return a list
#' @export
#'
#' @examples
#' source <- list(
#'   total   = 1,
#'   label   = c(2,3),
#'   unlabel = c(4,5))
#' pattern <- c("total", "total", "label", "total", "unlabel")
#' multiplyList(source, pattern)
#' # $total
#' # [1] 1
#' # 
#' # $total
#' # [1] 1
#' # 
#' # $label
#' # [1] 2 3
#' # 
#' # $total
#' # [1] 1
#' # 
#' # $unlabel
#' # [1] 4 5
multiplyList <- function(source, pattern) {
  res <- list()
  for (i in seq_along(pattern)) {
    res[[i]] <- source[[as.character(pattern[i])]]
  }
  names(res) <- pattern
  res
}


#' Leave only one item per name 
#' 
#' If the same name is used for several items (e.g. "total"),
#' only first item is left in the result.
#'
#' @param list a list to process
#' 
#' @return a shorter list
#' @export
#'
#' @examples
#' l <- list(
#'   total   = 1,
#'   total   = 1,
#'   unlabel = c(4,5))
#' shrinkList(l)
#' # $total
#' # [1] 1
#' # 
#' # $unlabel
#' # [1] 4 5
#' 
shrinkList <- function(list){
 list[unique(names(list))] 
}

assignList <- function(l, x){
 utils::relist(rep(x, length(unlist(l))), l) 
}

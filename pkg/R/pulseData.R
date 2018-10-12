
#' Create an object for pulse-change count data
#'
#' @param counts a matrix; column names correspond to sample names.
#' The columns in `counts` correspond to the rows in `conditions` argument.
#' @param conditions a data.frame;
#'   the first column corresponds to the conditions given in \code{formulas}.
#'   The order of rows corresponds to the columns (samples) in the
#'   `counts` argument.
#' @param formulas a list, created by \code{\link{MeanFormulas}}
#' @param formulaIndexes a list of lists (or of vectors);
#' defines indexes of formulas used for calculation of the expected read number.
#'
#' @param spikeins NULL (default) or a list of two items:
#'   - refGroup, a character, which defines the group which should be
#'     treated as a reference for normalisation
#'   - spikeLists, a list of character vectors with the spike-ins names and
#'     the same structure as `formulaIndexes`.
#' @param groups NULL (default) or a vector or a formula,
#'  e.g. ~ fraction + time.
#' If the normalisation factors must be recovered during fitting,
#' `groups` define the sets of the samples, which share the same normalisation
#' factors. Hence `groups` is relevant only, if there were no spike-ins provided.
#' In this case, one may assume that for the samples of the same nature, i.e.
#' same fraction and labelling time, the efficiency of the purification is the
#' same, which reduces the number of parameters to fit. However, it is possible
#' to treat every sample as an individual group, hence there will be no shared
#' parameters.
#'
#' If it is NULL, `groups` are derived from the first column of the `conditions`.
#' If a vector is provided, its elements correspond to the rows (samples) in the
#' `conditions`.
#' If, for example, there are 3 pull-down (2hr) samples purified with the protocol
#' "A", and 3 pull-down (2hr) samples from the protocol "B", one may assume different
#' efficiency of this protocols and reflect it in the `groups` argument by
#' introducing additional column `protocol` in the condition matrix, which results in
#' `groups = ~ fraction + time + protocol`.
#' Alternatively, one may manually create a vector like
#' `c("pull_down.2hr.A", "pull_down.2hr.B", ...)` with the order, corresponding
#' to the sample order in the `conditions`.
#'
#' @return an object of class "PulseData"
#'
#' It is a list with the following slots:
#'
#' - `user_conditions`, `user_formulas`, `counts` are the values of arguments
#' `conditions`, `formulas` and `counts`, provided to the call of `PulseData`
#' - `rawFormulas` is a list of initial formulas, evaluated at the
#' corresponding conditions (e.g. `time` in formulas is substituted with
#' its values in the `conditions$time`).
#' - `formulas` is a list of the compiled `rawFormulas`
#' - `formulaIndexes` is a list of integers (or vectors) with indexes of formulas
#'   used in estimation of the expression level in a given sample.
#'   The order of list items corresponds to the order of the samples in the
#'   `conditions` data.frame. See also `\link{addKnownToFormulas}`.
#' - `groups` is a vector with the names of the sample groups, which is used to
#'   calculate normalisation factors.
#' - `depthNormalisation` is a list of normalisation factors of the same structure as
#'   `formulaIndexes`. If no spike-ins are used, these values correspond
#'   to sequencing depth within a given group of samples according to the
#'   `groups` vector. For example, depth normalisation for a group "pull_down.2hr"
#'   of the pull-down samples after 2 hr of labelling. The relation between
#'   different groups, i.e. "total_fraction", "pull_down.2hr" etc., is
#'   not known and must be recovered during fitting as `normFactors` values.
#'   If spike-ins are provided, the relation between different fractions
#'   is recovered during the initialisation of the `PulseData` object and
#'   the values are written to the `depthNormalisation` slot.
#' - `interSampleCoeffs` is a list, which structure is used as a sekeleton for the
#'   normalisation factors, if no spike-ins were provided. For every group in
#'   `groups`, there is a corresponding list item (a number of a numeric vector).
#' - `interSampleIndexes` describes which normalisation factor to use during
#'   calculation of 
#'
#' @details
#' The `conditions` argument may include additional  columns, which
#' provide values for known parameters, such as time. Their name must be the
#' same as defined in formulas. For example, if a formula is defined as
#' `mu * exp(-d * time)` where `time` is the time point of the experiment,
#' the condition data.frame must contain a column named `time`, otherwise `time`
#' is treated as a parameter to fit!
#' @import methods
#' @export
#'
#' @examples
#'
#'
#'
#' formulaIndexes <- list(
#'   total_fraction = 'total',
#'   flow_through   = c('unlabelled', 'labelled'),
#'   pull_down      = c('labelled', 'unlabelled'))
#'
#' # Spike-ins definition for object creation
#' refGroup <- "total_fraction"
#'
#' labelled <- c("spike1", "spike2")
#' unlabelled <- c("spike3", "spike4")
#'
#' spikeLists <- list(
#' # total samples are normalised using all spike-ins
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
  e <- list()
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
    # if no groups are provided every sample is treated individually
    # this leads to deseq depth normalisation factors of 1
    # and estimation of normalisation factors by MLE
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
 print(table(x$conditions[,1]))
}

#' Calculate normalisation factors for columns in a matrix
#'
#' @param count_data a matrix; columns correspond to samples.
#'
#' @return a vector of doubles with normalisation factors ordered as columns in 
#'   the \code{count_data}
#'
#' @importFrom  stats median
#' @keywords  internal
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
    print(utils::head(x))
    stop("Can't normalise accross a condition.
         Too many zero expressed genes. ")
  }
  res
}

#' Estimate DESeq normalisation factors using spike-ins counts
#'
#' @param pd a \code{\link{PulseData}} object
#' @param refGroup a character name of a sample groups to use as a reference.
#' The name must be present in the spikeLists names list
#' @param spikeLists a named list with the same structure as formulaIndexes
#'
#' @return a list of lists of the normalisation coefficients for every sample
#' @keywords  internal
#'
normaliseWithSpikeIns <- function(pd, refGroup, spikeLists){
  refSpikes <- unlist(spikeLists[[refGroup]])
  # Subset only spike-ins contained in the reference 
  spikeCounts <- pd$counts[refSpikes,]
  spikeLists <- lapply(spikeLists,
                       function(l) {
                         lapply(l, function(spikes)
                           spikes[spikes %in% refSpikes])
                       })
  # create DESeq reference virtual sample (geo-means of counts) 
  # and compute DESeq factors using spike-ins common with the reference
  superSample <- rowMeans(
    log(pd$counts[refSpikes, pd$conditions[,1] == refGroup]))
  lapply(seq_along(pd$conditions[, 1]),
         function(i) {
           sampleSpikes <- spikeLists[[as.character(pd$conditions[i, 1])]]
           vapply(sampleSpikes, function(spikes) {
             deseq(spikeCounts[spikes, i], superSample[spikes])
           }, double(1))
         })
}


#' Performs sequencing depth normalisation using the DESeq procedure.
#'
#' @param pd a \code{\link{PulseData}} object 
#' @param groups a vector for splitting objects to groups
#'
#' @return a vector with the coefficient for every sample
#' @keywords  internal
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
#' @keywords  internal
#' 
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
#'  for expected levels in a given fraction
#'  (e.g. describing amounts of "labelled", "total" RNA etc.)
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
#'
addKnownToFormulas <- function(formulas, formulaIndexes, conditions) {
  uc <- unique(conditions)
  uc <- uc[order(uc[,1]),,drop = FALSE]
  newIndexes <- list()
  newForms <- list()
  for (i in order(uc[, 1])) {
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
#' @keywords  internal
#'
makeGroups <- function(pd, normGroups) {
  # generate a list of normalisation coefficients with a proper structure
  # names are ordered according to the first columns of conditions
  groupNames <- unique(normGroups[order(pd$conditions[,1])])
  normCoeffs <- pd$formulaIndexes[match(groupNames, normGroups)]
  names(normCoeffs) <- groupNames
  normCoeffs[order(normGroups)]
  # all the normalisation coefficients are numbered according 
  # to their appearance in the flatten list `unlist(normCoeffs)`
  normCoeffs <- utils::relist(seq_along(unlist(normCoeffs)), normCoeffs)
  # normCoeffIndexes contains items for every sample
  normCoeffIndexes <- multiplyList(normCoeffs, normGroups)
  normCoeffs <- assignList(normCoeffs, 1)
  list(normCoeffs       = normCoeffs,
       normCoeffIndexes = normCoeffIndexes)
}

# Returns a list of the same structure as nameList
# with the names substituted by their index in the nameVector
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
  evaled <- do.call(cbind, lapply(formulas, eval, env = par))
  # if normFactors are specified only once per condition --> multiply
  if (length(normFactors) == length(unique(conditions)[,1])) {
    normFactors <- multiplyList(normFactors, conditions[,1])
  }
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
shrinkList <- function(list){
 list[unique(names(list))] 
}

assignList <- function(l, x){
 utils::relist(rep(x, length(unlist(l))), l) 
}

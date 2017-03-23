#' An example of count data for an experiment without spike-ins
#' 
#' @details 
#' The data set contains simulation of an experiment which measures three 
#' fractions, namely, coded as 'A_fraction', 'B_fraction' and 'C_fraction'.
#' There are three types of quantities, 'A', 'B' and 'C' for which a kinetic
#' model is defined, see `formulas` element in the data set:
#' 
#' (A = a * p, B =  a * b^time, C = a * (1 - b^time)),
#' 
#' where a, b are gene-specific parameters which are unknown, p is 
#' a vector of known gene-specific parameters (e.g. probability of capture
#' depending on the uridine content etc.)
#' 
#' The data are generated for 3 replicates, 3 different time points and for
#' 10 different genes, see elements `counts` and `conditions`.
#' 
#' The model considers possibility of cross-contamination with different types
#' of RNA, which is described by `formulaIndexes`. In this case,
#' the mean read number for a gene is a linear combination of the described
#' RNA types with weights defined in `normFactors`.
#' 
#' Due to different amounts of RNA in different conditions, the normalisation
#' factors are defined for groups of the samples, hence they are shared values
#' for a given sample group. The grouping is defined by `groups` element of the
#' data set.
#' 
#' The true parameter values, which are used for data simulation are in 
#' `par` element. This also includes the size parameter for the negative
#' binomial distribution.
#' 
#' 
#' @format A list containing simulated data, model and parameter information.
#' 
#' - `formulas` describe the model for mean read number
#' - `counts` contain the simulated data
#' - `conditions` is a data.frame describing sample time point and type    
#' - `fractions` is a character vector to divide samples into groups for 
#'   normalisation
#' - `formulaIndexes` is a list of formula names (see `formulas`), which are
#'   used in linear combination with coefficients in `normFactors` to calculate
#'   mean read number
#' - `normFactors` are coefficients which relate quantities of RNA between 
#'   different groups and can be used as normalisation coefficients for 
#'   calculation of mean read numbers using `formulas` and `formulaIndexes`   
#' - `par` is a list of model parameters 
#'
#'
"pulseRFractionData"



#' An example of count data for an experiment with spike-ins
#'
#' @details 
#' The data set contains simulation of an experiment which measures three 
#' fractions, namely, coded as 'A_fraction', 'B_fraction' and 'C_fraction'.
#' There are three types of quantities, 'A', 'B' and 'C' for which a kinetic
#' model is defined, see `formulas` element in the data set:
#' 
#'   A = a,
#'   B =  a * b ^ time,
#'   C = alpha * a * (1 - b ^ time)
#' 
#' where a, b are gene-specific parameters which are unknown, 
#' alpha is a parameter which is shared between all genes.
#' 
#' The data are generated for 3 replicates, 3 different time points and for
#' 10 different genes, see elements `counts` and `conditions`.
#' 
#' The model considers possibility of cross-contamination with different types
#' of RNA, which is described by `formulaIndexes` simply as
#' 
#' formulaIndexes <- list(
#'   A_fraction = 'A',
#'   B_fraction = c('B', 'C'),
#'   C_fraction = c('B', 'C'))
#'   
#' In this case,
#' the mean read number for a gene is a linear combination of the described
#' RNA types with weights defined in `normFactors`.
#' 
#' Spike-ins counts are generated in order to recover normalistion coefficients
#' which describe how read counts in different samples relate to each other.
#' Different spike-ins correspond to different types of RNA (e.g. labelled and 
#' unlabelled) and the rule for this  relations are defined in the 
#' the `spikeins` element for this data set.
#' 
#' The true normalisation coefficients which were used for data simulation
#' are contained in the `allNormFactors` list.
#' 
#' The true parameter values, which are used for data simulation are in 
#' `par` element. This also includes the size parameter for the negative
#' binomial distribution.
#' 
#' @format A list containing simulated data, model and parameter information.
#' 
#' - `formulas` describe the model for mean read number
#' - `counts` contain the simulated data
#' - `conditions` is a data.frame describing sample time point and type    
#' - `fractions` is a character vector to divide samples into groups for 
#'   normalisation
#' - `formulaIndexes` is a list of formula names (see `formulas`), which are
#'   used in linear combination with coefficients in `normFactors` to calculate
#'   mean read number
#' - `spikeins` is a list with two elements:
#'   - `refGroup` is the name of the sample which to use as a reference sample
#'   - `spikeLists` is a list of names of spike-ins which are to be used
#'     for normalisation. It has the same structure as `formulaIndexes` 
#' - `allNormFactors` is a list of true coefficients which were used for
#'   spike-ins and sample counts normalisations before simulation
#' - `par` is a list of model parameters 
#' 
#'
"pulseRSpikeinsData"
#' An example of count data for an experiment without spike-ins
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
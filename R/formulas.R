
toLanguage <- function(x) {
  if (!is.language(x))
    parse(text = x)[[1]]
  else
    x
}

#' Create a formula for RNA degradation
#' 
#' All the arguments must be characters or language objects 
#' (expression, call or name).
#' 
#' @param x  initial concentration.
#' @param d  a degradation rate.
#' @param t  a longitude of the modelled period.
#'
#' @return an expression for the calculation of the RNA level after
#' degradation during time \verb{t}.
#' @export
#'
#' @examples
#' x <- quote(mu_0)
#' mu <- quote(mu_new)
#' d <- "degradation_rate"
#' t <- "t_labelling"
#' degrade(x, d, t)
#' # mu_0 * exp(-degradation_rate * t_labelling)
#' 
degrade <- function(x, d, t){
  args <- as.list(match.call()[-1])
  args <- lapply(args, eval, env=parent.frame())
  args <- lapply(args, toLanguage)
  do.call(degrade_, args) 
}

#' Create a formula for RNA degradation
#' 
#' Implements a non-standard evaluation version.
#'
#' @inheritParams degrade
#'
#' @inheritSection return 
#' @export
#'
#' @examples
#' degrade_(a,b,c)
#' # a * exp(-b * c)
#' 
degrade_ <- function(x, d, t) {
  args <- as.list(match.call()[-1])
  bquote(.(x) * exp(-.(d) * .(t)), args)
}

#' Creates a formula which describe evolution of RNA concentration
#'
#' All the arguments must be characters or language objects 
#' (expression, call or name).
#' 
#' @param x  initial concentration.
#' @param mu new steady-state level.
#' @param d  a degradation rate.
#' @param t  a longitude of the modelled period.
#'
#'
#' @return an expression for the calculation of the RNA level after time
#' \verb{t}.
#' @export
#'
#' @examples
#' x <- quote(mu_0)
#' mu <- quote(mu_new)
#' d <- "degradation_rate"
#' t <- "t_labelling"
#' grow(x, mu, d, t)
#' # mu_new - (mu_new - mu_0) * exp(-degradation_rate * t_labelling)
grow <- function(x, mu, d, t) {
  args <- as.list(match.call()[-1])
  args <- lapply(args, eval, env=parent.frame())
  args <- lapply(args, toLanguage)
  do.call(grow_, args) 
}


#' Creates a formula which describe evolution of RNA concentration
#' 
#' This implement the non-standard evaluation version of the 
#' \code{\link{grow}} function.
#'
#' @inheritParams grow
#' @return an expression for the calculation of the RNA level.
#' @export
#'
#' @examples
#' grow_(a,b,c,d)
#' # b - (b - a) * exp(-c * d)
#' 
grow_ <- function(x, mu, d, t){
  args <- as.list(match.call()[-1])
  bquote(.(mu) - (.(mu) - .(x)) * exp(-.(d) * .(t)), args)
}


`%*%` <- function(lhs, rhs){
  parent <- parent.frame()
  env <- new.env(parent = parent)  
}

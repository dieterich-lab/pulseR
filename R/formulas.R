
toLanguage <- function(x) {
  if (!is.language(x))
    parse(text = x)[[1]]
  else
    x
}

#' Create a formula object for the initial RNA level.
#'
#' @param x a character or a language object (an expression, a call
#'        or a name)
#'
#' @return an expression for the initial RNA level
#' @export
#'
#' @examples
#' x <- amount("mu")
#' degrade_(x, d_rate, time)
#' # mu_0 * exp(-d_rate * time)
amount <- function(x){
 x <- toLanguage(x)
 x
}

#'  Create a formula object for the initial RNA level.
#'
#'  A non-standard evaluation version of the \code{\link{amount}} 
#'  function.
#'  
#' @param x an expression level
#'
#' @return an expression for the initial RNA level
#' @export
#'
amount_ <- function(x){
 substitute(x) 
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
#' x <- amount("mu_0")
#' mu <- amount_(mu_new)
#' d <- "degradation_rate"
#' t <- "t_labelling"
#' degrade(x, d, t)
#' # mu_0 * exp(-degradation_rate * t_labelling)
#' 
degrade <- function(x, d, t){
  args <- as.list(match.call()[-1])
  args <- lapply(args, eval, env = parent.frame())
  args <- lapply(args, toLanguage)
  args$x <- quote(x)
  do.call(degrade_, args) 
}

#' Create a formula for RNA degradation
#' 
#' Implements a non-standard evaluation version.
#'
#' @inheritParams degrade
#'
#' @return an expression for the calculation of the RNA level after
#' degradation during time \verb{t}.
#' @export
#'
#' @examples
#' x <- amount("a")
#' degrade_(x,b,c)
#' # a * exp(-b * c)
#' 
degrade_ <- function(x, d, t) {
  args <- as.list(match.call()[-1])
  args$x <- x
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
#' x <- amount_(mu_0)
#' mu <- amount("mu_new")
#' d <- "degradation_rate"
#' t <- "t_labelling"
#' grow(x, mu, d, t)
#' # mu_new - (mu_new - mu_0) * exp(-degradation_rate * t_labelling)
grow <- function(x, mu, d, t) {
  args <- as.list(match.call()[-1])
  args <- lapply(args, eval, env=parent.frame())
  args <- lapply(args, toLanguage)
  args$x <- quote(x)
  do.call(grow_, args) 
}

grow0 <- function(mu, d, t) {
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
#' x <- amount("a")
#' grow_(x,b,c,d)
#' # b - (b - a) * exp(-c * d)
#' 
grow_ <- function(x, mu, d, t){
  args <- as.list(match.call()[-1])
  if(!missing(x)){
		args$x <- x
		bquote(.(mu) - (.(mu) - .(x)) * exp(-.(d) * .(t)), args)
  } else {
		bquote(.(mu) - .(mu) * exp(-.(d) * .(t)), args)
  }
}


`%$%` <- function(lhs, rhs){
  parent <- parent.frame()
  env <- new.env(parent = parent)  
}

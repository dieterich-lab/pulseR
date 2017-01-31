
degrade_ <- function(x, d, t){
	substitute( x * exp(-d*t), 
		list(x=as.name(x),
				 d=as.name(d),
				 t=as.name(t)))
}

degrade <- function(x, d, t){
	substitute( x * exp(-d*t))
}

grow_ <- function(x, mu, d, t){
	substitute(mu - (mu - x) * exp(-d*t), 
		list(x=as.name(x),
				 d=as.name(d),
				 t=as.name(t)))
}
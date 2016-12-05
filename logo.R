library(ggplot2)
library(magrittr)


#q0 <- ggplot() +
#  geom_polygon(aes(x=c(-R,R,R,-R), y=c(R,R,-R,-R)), fill="black", alpha=.8) +
#  theme_void()
#q0 <- ggplot() + theme_void() +
#   theme(panel.background = element_rect(fill=alpha("blue",.8)))

f <- function(x,s,R, r=1) {
  x <- x/r
  res <- (exp(-(x+s/2)^2) - exp(-(x-s/2)^2)*.9)*R
  res *R /max(res)
}

DNA <- function(a,b,y,size, n, phi0,s,decay){
  x <- seq(from=a, to=b, length.out=1000)
  A <- rep(size, length(x))
  o <- round(length(x) * decay) :length(x)
  A[o] <- size* sqrt(b-x[o])/sqrt(b-x[o][1])
  p <- 2*pi*n/abs(b-a)*x + phi0
  y1 <- A*  sin(p) + y
  p2 <- p + s/abs(a-b)*2*pi*n
  y2 <- A*  sin(p2) + y
  y1[p < p2[1]] <- NA
  data.frame(x=x,y1=y1,y2=y2)
}


## Ploting
plotCircle <- function(color="black", alpha=.8){
  phi <- 2*pi * 1:101/100
  circle <- data.frame(x=R*sin(phi), y=R*cos(phi))
  q0 <- ggplot() + theme_void()
  if(missing(color))
   q0 <- q0 +
    geom_polygon(data=circle, aes(x=x,y=y), fill=color, alpha=alpha ) +
   theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(fill = "transparent",colour = NA)
  ) else
   q0 <- q0 +
    geom_polygon(data=circle, aes(x=x,y=y), fill=color, alpha=0) +
   theme(
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(fill = alpha(color, alpha),colour = NA)
    )
  q0
}
plotCircle("black")

plotDNA <- function(q0, t=1, color="green"){
  n <- dim(DNA1)[1]
  nas <- is.na(DNA1$y1)
  q1 <- q0 +
    geom_line(data=DNA1[1:round(t*n),],
              aes(x=x, y=y2), color=color, alpha=.7, size=lineWidth) +
    geom_line(data=DNA1[!nas & 1:n < round(t*n),],
              aes(x=x, y=y1), color=color, alpha=.8, size=lineWidth)
  if(t < 1){
    p <- DNA1[round(t*n),]
    q1 <- AddPoint(q1, p$x, p$y1)
    q1 <- AddPoint(q1, p$x, p$y2)
  }
  q1
}

plotPulse <- function(q1, t=1, color="green"){
  n <- dim(pulse)[1]
  q2 <- q1 + geom_line(data=pulse[1:round(t*n),], aes(x=x,y=y),
                      color=color, alpha=.9, size=lineWidth)
  if(t < 1){
    p <- pulse[round(t*n),]
    q2 <- AddPoint(q2, p$x, p$y)
  }
  q2
}

AddText <- function(q, color="white"){
  q3 <-  q + geom_text(aes(x=4.3, y=1.0, label="pulseR"),
                        color=color, size=28, family="Ubuntu Mono"
                        , fontface="italic"
                      )
  q3
}

AddPoint <- function(q,x,y){
  d <- data.frame(x=x,y=y)
  q +
  geom_point(data=d, aes(x=x,y=y),
                 color="green",
                 fill="green",
                 size=15,
                 shape=21,
                 stroke=0,
                 alpha=.2) +
  geom_point(data=d, aes(x=x,y=y),
                 color="green",
                 size=5) +
  geom_point(data=d, aes(x=x,y=y),
                 color="white",
                 size=1)
}

lineWidth <- 2
R <- 10
shift <- -3.4
phi0 <- -2.1

shiftPulse <- -1.7
pulse <- seq(shiftPulse,R - shift,.1) %>% data.frame(x=.,y=f(.,.5,R*.9,.8))
pulse$x <- pulse$x + shift - shiftPulse
pulse <- pulse[ pulse$x < R,]
pulse$y <- pulse$y - pulse$y[1]


DNA1 <- DNA(a=-6,b=0, y=0, size=3,n=2,phi0=phi0,s=.7,decay=.2)
DNA1$x <- DNA1$x + shift

## plain logo
q <- plotCircle()
q <- AddText(q)
q <- plotDNA(q)
q <- plotPulse(q)
ggsave(filename="logo.png", plot=q, height=7,width=7)

q <- plotCircle("white")
q <- AddText(q, "black")
q <- plotDNA(q, color="darkgreen")
q <- plotPulse(q, color="darkgreen")
ggsave(filename="logo_white.png", plot=q, height=7,width=7)

## animation

trace <- function(bgcolor, dnacolor="green",textcolor="white", alpha=.8){
  if(missing(bgcolor))
    q <- plotCircle()
  else
    q <- plotCircle(bgcolor, alpha)
  q <- AddText(q, textcolor)
  t <- seq(from=0, to=1, by=1/35)
  lapply(t, function(x){
     print(plotDNA(q,t=x, color=dnacolor))})
  q <- plotDNA(q, color=dnacolor)

  frac <- .3
  t <- seq(from=0, to=frac, by=1/70)
  lapply(t, function(x){
     print(plotPulse(q,t=x, color=dnacolor))})
  t <- seq(from=frac, to=1, by=1/70)
  lapply(t, function(x){
     print(plotPulse(q,t=x, color=dnacolor))})
}

saveGIF(trace(), interval=.02, movie.name="logo.gif")
saveGIF(trace("black"), interval=.02,
        movie.name="logo_background.gif")
saveGIF(trace("white", "darkgreen", "black"), interval=.02,
        movie.name="logo__white_background.gif")
library(ggplot2)
library(magrittr)

lineWidth <- 2

phi <- 2*pi * 1:101/100
R <- 10
circle <- data.frame(x=R*sin(phi), y=R*cos(phi))
q0 <- ggplot() +  geom_polygon(data=circle, aes(x=x,y=y)) + theme_void()

#q0 <- ggplot() +
#  geom_polygon(aes(x=c(-R,R,R,-R), y=c(R,R,-R,-R)), fill="black", alpha=.8) +
#  theme_void()
#q0 <- ggplot() + theme_void() +
#   theme(panel.background = element_rect(fill=alpha("blue",.8)))

f <- function(x,s,R, r=1) {
  x <- x/r
  (exp(-(x+s/2)^2) - exp(-(x-s/2)^2))*R
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

shift <- -3.4
phi0 <- -2.1

shiftPulse <- -1.7
pulse <- seq(shiftPulse,R - shift,.1) %>% data.frame(x=.,y=f(.,.5,2*R,.8))
pulse$x <- pulse$x + shift - shiftPulse
pulse <- pulse[ pulse$x < R,]


DNA1 <- DNA(a=-6,b=0, y=pulse$y[1]+.1,size=3,n=2,phi0=phi0,s=.7,decay=.2)
DNA1$x <- DNA1$x + shift

q1 <- q0 +
  geom_line(data=DNA1, aes(x=x, y=y1),color="green", alpha=.7, size=lineWidth) +
  geom_line(data=DNA1, aes(x=x, y=y2), color="green", alpha=.8, size=lineWidth)

q2 <- q1 + geom_line(data=pulse, aes(x=x,y=y),
                     color="green", alpha=.9, size=lineWidth)

q3 <-  q2 + geom_text(aes(x=4.3, y=1.5, label="pulseR"),
                      color="white", size=25, family="NimbusSan"
                      , fontface="italic"
                     )
q3
ggsave(filename="logo.png", plot=q3, height=7,width=7)
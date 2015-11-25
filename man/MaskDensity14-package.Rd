\name{MaskDensity14-package}
\alias{MaskDensity14-package}
\alias{MaskDensity14}
\docType{package}
\title{
An R package for density approximant based on noise multiplied data
}
\description{
This package is used to estimate the density approximant of the original data,
in the meanwhile  to produce synthetic data of the original data,
 based on noise multiplied data.

}
\details{
\tabular{ll}{
Package: \tab MaskDensity14\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-20\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Mark and Yan-Xia

Maintainer: Yan-Xia Lin <yanxia@uow.edu.au>
~~ The author and/or maintainer of the package ~~
}
\references{
no reference
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ MaskDensity14 }


\examples{
# Example for continuous data

# Data provider would have the real data, in this case we are just simulating some continuous data: 
set.seed(123)
n=10000
rmulti<-function(n,mean,sd,p)
   {
     x<-rnorm(n,0,1)
     k<- length(mean)
     u<-sample(1:k, size=n, prob=p, replace=TRUE)
     for(i in 1:k)
         x[u==i]<-mean[i]+sd[i]*x[u==i]
         
     return(x)
   }
   y<-rmulti(n=10000, mean=c(30,50),sd=c(4,2), p=c(0.3,0.7))
       # y is a sample drawn from Y
   noise <- rmulti(n=10000, mean=c(80,100), sd=c(5,3), p=c(0.6,0.4))
       # noise is a sample drawn from C

# Data provider now masks the real data created above:
   library(MaskDensity14)
   ymask<- mask(y, noisefile="noise.bin", a1=min(y), b1=max(y), noise)
   write(ymask$ystar, "ystar.dat")
   
# Data user unmasks the masked data:   
   library(MaskDensity14)
   ystar<-scan("ystar.dat")
   y1<-unmask(ystar, noisefile="noise.bin")
   plot(density(y1$unmaskedVariable), main="density(ymask)", xlab="y")
   
# Example for categorical data

# Data provider would have the real data, in this case we are just simulating some categorical data: 
set.seed(124)
n<-2000
a<-170
b<-80
y<-rbinom(n, 1, 0.1)+1
noise<-(a+b)/2+ sqrt(1+(a-b)^2/4)*rnorm(n, 0,1)
noise<-abs(noise)

# Data provider now masks the real data created above:
library(MaskDensity14)
ymask<-mask(factor(y), noisefile="noise.bin", noise, a1=0,b1=3)
      # using factor(y) because y is a categorical variable
write(ymask$ystar, "ystar.dat")

# Data user unmasks the masked data:   
library(MaskDensity14)
ystar<-scan("ystar.dat")
y1 <- unmask(ystar, noisefile="noise.bin")
unmaskY<-y1$unmaskedVariable  # synthetic data
mass_function<-y1$prob  # estimated mass function



}

\name{mask}
\alias{mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
masking data
}
\description{
mask the original data by noise
}
\usage{
mask(y, noisefile, noise = createNoise(length(y)), a1 = min(y), b1 = max(y), maxorder = 100, EPS = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the original data
}
  \item{noisefile}{
a file containing the information of noise
}
  \item{noise}{
the noise used to mask data
}
  \item{a1}{
the a_{basic} mentioned in paper
}
  \item{b1}{
the b_{basic} mentioned in paper
}
  \item{maxorder}{
the maximum order of moments used in the density approximant
}
  \item{EPS}{
the criterion used in order selection
}
}
\details{
no details needed
}
\value{

 \item{ystar }{masked data}
 \item{noisefile}{a file containing the information of noise}

%% ...
}
\references{
no reference
}
\author{
Mark and Yan-Xia
}
\note{
no further notes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~



\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

y<-rnorm(100,0,1)
ymask<-mask(y, noisefile="noise.bin", noise=createNoise(length(y)), a1=min(y), 
     b1=max(y), maxorder=100, EPS=1e-06)
write(ymask$ystar, "ystar.dat")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mask_data  }
 
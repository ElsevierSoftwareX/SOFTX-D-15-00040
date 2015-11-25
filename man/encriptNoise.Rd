\name{encriptNoise}
\alias{encriptNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
encript noise
}
\description{
to create a binary file for noise
}
\usage{
encriptNoise(noise, a, b, maxorder, levels, EPS, noisefile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noise}{
noise used to mask data
}
  \item{a}{
lower boundary
}
  \item{b}{
upper boundary
}
  \item{maxorder}{
determine the maximum order of the moment in Provost
}
  \item{levels}{
levels for categorical variable
}
  \item{EPS}{
criterion used to stop the number of moments in provost
}
  \item{noisefile}{
a binary noise file
}
}
\details{
no details required
}
\value{
binary code

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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

c2<-rnorm(100,0,1)
c2<-abs(c2)
a<-0
b<-2
maxorder<-10
lvls<-NULL
EPS<-1e-06
encriptNoise(c2,a,b,maxorder,lvls,EPS,"noisefile")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{noisefile }
 
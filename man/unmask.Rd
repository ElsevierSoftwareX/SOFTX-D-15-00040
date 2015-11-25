\name{unmask}
\alias{unmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the main function
}
\description{
The final code for the density approximant of the original data
}
\usage{
unmask(ystar, noisefile, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ystar}{
masked data
}
  \item{noisefile}{
a file containing the information of noise
}
  \item{noise}{
noise used to mask the data y
}
}
\details{
no details needed
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{unmaskedVariable}{synthetic data for y}
  \item{prob}{estimated mass function for y}
  \item{meanOfNoise}{the sample mean of noise}
  \item{meanOfSquaredNoise}{the sample mean of the second moment of the noise}
%% ...
}
\references{
no reference needed
}
\author{
Yan-Xia
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
ymask<-mask(y,noisefile="noise.bin", a1=min(y), b1=max(y), noise=createNoise(length(y)),
       maxorder=100, EPS=1e-06)
write(ymask$ystar, "ystar.dat")
ystar<-scan("ystar.dat")
y1<-unmask(ystar,noisefile="noise.bin", noise=noise)
print(y1$unmaskedVariable)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{density_approximant  }


\name{createNoise}
\alias{createNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a noise
}
\description{
This function is used to create a noise with 5-modal mixture normal distribution
}
\usage{
createNoise(n, mean = rep(NA, 5), sd = rep(1, length(mean)), prob = rep(1, length(mean))/length(mean))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the size of noise
}
  \item{mean}{
the size of noise
}
  \item{sd}{
a vector of the standard error of the 5-modal noise
}
  \item{prob}{
a vector of the propotion of the sample from the 5-modal distribution
}
}
\details{
no detalis required
}
\value{
a real vector

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
no references
}
\author{
Mark and Yan-Xia
}
\note{
no further notes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

y<-rnorm(100)
noise<- createNoise(length(y))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ noise }


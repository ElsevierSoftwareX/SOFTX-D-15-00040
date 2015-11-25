\name{moments}
\alias{moments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
moment
}
\description{
calculate the values of the moments of Y
}
\usage{
moments(y, order = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a sample from Y
}
  \item{order}{
the maximum of order of the moments used
}
}
\details{
no details needed
}
\value{
a vector

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
no reference needed
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
print(moments(y, order=8))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{moments  }
 

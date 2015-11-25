\name{calc_muX}
\alias{calc_muX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the moments
}
\description{
Calculate mu_X(j) in equation (9) in Provost paper
}
\usage{
calc_muX(muY, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{muY}{
moment of Y
}
  \item{a}{
lower boundary
}
  \item{b}{
upper boundary
}
}
\details{
there are no more  details required
}
\value{
A vector
}
\references{
Moment-Based Density Approximants by S.B. Provost
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

a<-2
b<-5
muY<-c(3,3.5,5)
print(calc_muX(muY, a,b))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ moment }

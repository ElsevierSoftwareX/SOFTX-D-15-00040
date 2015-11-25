\name{fY_Rmask}
\alias{fY_Rmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
density of Y
}
\description{
evaluate the density function for the general case
}
\usage{
fY_Rmask(y, muY, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
values of Y
}
  \item{muY}{
a vector with entries of moments of Y
}
  \item{a}{
lower boundary of Y
}
  \item{b}{
upper boundary of Y
}
}
\details{
no details needed
}
\value{
a real number

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Equation (12) in Provost's paper
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

y<-4
muY<-c(3,3,5,5)
a<-2
b<-5
print(fY_Rmask(y,muY,a,b))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{lambda_Rmask}
\alias{lambda_Rmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
expectation of Legendre polynomial
}
\description{
calculate the lambda function given by Equation (3) in Provost's paper
}

\usage{
lambda_Rmask(muX, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{muX}{
a vector of the moments of X
}
  \item{k}{
the upper order of the moments used in the Legendre polynomial
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
Equation (3) in provost's paper
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

muX<-c(1,2,3,4,5)
k<-4
print(lambda_Rmask(muX,k))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{lambda }
 
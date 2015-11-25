\name{findOrder_Rmask}
\alias{findOrder_Rmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
order determination
}
\description{
Determine the upper order of moment in the Legendre polynomial function
}
\usage{
findOrder_Rmask(ystar, noise, a, b, maxorder = 100, EPS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ystar}{
masked data
}
  \item{noise}{
noise information
}
  \item{a}{
lower boundary
}
  \item{b}{
upper boundary
}
  \item{maxorder}{
the maximum order of the moments used in the analysis
}
  \item{EPS}{
criterion for stopping
}
}
\details{
no details needed
}
\value{

integer

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
the idea is described in the paper
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

## The function is currently defined as
y<-rnorm(100,0,1)
noise<-abs(rnorm(100,2,4))
ystar<- noise*y
a<- -2
b<- 2
order<-findOrder_Rmask(ystar,noise,a,b,maxorder=100, EPS=1e-06)
print(order)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

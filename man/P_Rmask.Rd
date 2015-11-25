\name{P_Rmask}
\alias{P_Rmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Legendre polynomial
}
\description{
An R code for a Legendre polynormal of degree k in x
}
\usage{
P_Rmask(x, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
variable in the polynomial
}
  \item{k}{
the order in the polynormal
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
Equation (2) in Provost paper
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

x<-5
k<-4
print(P_Rmask(x,k))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Legendre-polynomial  }


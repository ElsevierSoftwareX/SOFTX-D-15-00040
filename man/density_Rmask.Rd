\name{density_Rmask}
\alias{density_Rmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density function
}
\description{
estimate density function based on masked data
\usage{
density_Rmask(moments, a, b, n = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{moments}{
moments information
}
  \item{a}{
lower boundary
  \item{b}{
upper boundary
}
  \item{n}{
scale on x-axis
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

y<-rnorm(100, 0,1)
val <- 1
for(k in 1:8){
   val<- c(val, mean(y^k))
   }
a<- -2
b<- 2
Provost<- density_Rmask(val, a,b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{intcoef}
\alias{intcoef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Bilinear interpolation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intcoef(x1, x2, x1l, x1u, x2l, x2u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
}
  \item{x1l}{
%%     ~~Describe \code{x1l} here~~
}
  \item{x1u}{
%%     ~~Describe \code{x1u} here~~
}
  \item{x2l}{
%%     ~~Describe \code{x2l} here~~
}
  \item{x2u}{
%%     ~~Describe \code{x2u} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x1, x2, x1l, x1u, x2l, x2u) 
{
    X <- matrix(data = c(1, 0, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, 
        1, -1, -1, 1), nrow = 4, ncol = 4, byrow = TRUE)
    x11 <- (x1 - x1l)/(x1u - x1l)
    x22 <- (x2 - x2l)/(x2u - x2l)
    h <- matrix(data = c(1, x11, x22, x11 * x22), nrow = 1, ncol = 4) \%*\% 
        X
    h
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

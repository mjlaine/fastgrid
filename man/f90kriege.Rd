\name{f90kriege}
\alias{f90kriege}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
fortran interface to Krieging prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
f90kriege(x, y, b, grid, cy, cgrid, covpars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
}
  \item{cy}{
%%     ~~Describe \code{cy} here~~
}
  \item{cgrid}{
%%     ~~Describe \code{cgrid} here~~
}
  \item{covpars}{
%%     ~~Describe \code{covpars} here~~
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
function (x, y, b, grid, cy, cgrid, covpars) 
{
    nobs <- nrow(x)
    npar <- ncol(x)
    ngrid <- nrow(grid)
    if (nrow(y) != nobs | ncol(y) != 1 | nrow(cy) != nobs | ncol(grid) != 
        npar | nrow(cgrid) != ngrid | ncol(cy) != 2 | ncol(cgrid) != 
        2 | !is.numeric(covpars) | length(covpars) < 2) {
        stop("input dimensions do not match")
    }
    ypred <- .Fortran("kriegepred", x = as.double(x), y = as.double(y), 
        b = as.double(b), grid = as.double(grid), ypred = double(length = ngrid), 
        cy = as.double(cy), cgrid = as.double(cgrid), nobs = as.integer(nobs), 
        npar = as.integer(npar), ngrid = as.integer(ngrid), covpars = as.double(covpars))$ypred
    if (b[1, 1] == -1) {
        stop("Problems with the covariance information, no prediction done")
    }
    return(ypred)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

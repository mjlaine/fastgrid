\name{fastkriege}
\alias{fastkriege}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Faster Kriging
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates Kriging predictions
}
\usage{
fastkriege(trend_model, data, grid, cov.pars, bg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trend_model}{
%%     ~~Describe \code{trend_model} here~~
formula to build the regression matrices
}
  \item{data}{
%%     ~~Describe \code{data} here~~
data containing station temperatures as 
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
definition of model grid, must contain longitude, latitude coordinates and variables used in trend_model
}
  \item{cov.pars}{
%%     ~~Describe \code{cov.pars} here~~
cov.pars = c(sigmasq,phi,tausq)
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
Optional background field of same dimensions as grid.
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
Marko Laine
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{eagamma_cpp}
\alias{eagamma_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eagamma_cpp(nS, sS, tS, xS, yS, LS, US)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nS}{
%%     ~~Describe \code{nS} here~~
}
  \item{sS}{
%%     ~~Describe \code{sS} here~~
}
  \item{tS}{
%%     ~~Describe \code{tS} here~~
}
  \item{xS}{
%%     ~~Describe \code{xS} here~~
}
  \item{yS}{
%%     ~~Describe \code{yS} here~~
}
  \item{LS}{
%%     ~~Describe \code{LS} here~~
}
  \item{US}{
%%     ~~Describe \code{US} here~~
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
function (nS, sS, tS, xS, yS, LS, US) 
{
    .Call("_scale_eagamma_cpp", PACKAGE = "scale", nS, sS, tS, 
        xS, yS, LS, US)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

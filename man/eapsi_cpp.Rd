\name{eapsi_cpp}
\alias{eapsi_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eapsi_cpp(jS, sS, tS, mS, xoyS, uS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jS}{
%%     ~~Describe \code{jS} here~~
}
  \item{sS}{
%%     ~~Describe \code{sS} here~~
}
  \item{tS}{
%%     ~~Describe \code{tS} here~~
}
  \item{mS}{
%%     ~~Describe \code{mS} here~~
}
  \item{xoyS}{
%%     ~~Describe \code{xoyS} here~~
}
  \item{uS}{
%%     ~~Describe \code{uS} here~~
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
function (jS, sS, tS, mS, xoyS, uS) 
{
    .Call("_scale_eapsi_cpp", PACKAGE = "scale", jS, sS, tS, 
        mS, xoyS, uS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

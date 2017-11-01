\name{ac.phi}
\alias{ac.phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ac.phi(eta.pos, dapts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta.pos}{
%%     ~~Describe \code{eta.pos} here~~
}
  \item{dapts}{
%%     ~~Describe \code{dapts} here~~
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
function (eta.pos, dapts) 
{
    beta.pos <- un.transform(eta.pos)
    data.evaluation <- data.eval(beta.pos, 1:dsz)
    (data.evaluation$grad.log.pi \%*\% t(data.evaluation$grad.log.pi) + 
        sum(data.evaluation$lap.log.pi))/2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{datum}
\alias{datum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
datum(datum.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datum.idx}{
%%     ~~Describe \code{datum.idx} here~~
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
function (datum.idx) 
{
    old.seed <- .Random.seed
    on.exit({
        .Random.seed <<- old.seed
    })
    seed.indexing <- datum.idx\%/\%2147483647
    datum.seeding <- datum.idx\%\%2147483647
    set.seed(datum.seeding)
    runif((2 * dimen - 1) * seed.indexing)
    datum.x <- c(1, rtnorm(dimen - 1, mean = 0, sd = 1, lower = -design.thresh, 
        upper = design.thresh))
    datum.z <- exp(sum(beta.truth * datum.x))
    datum.y <- rbinom(1, 1, datum.z/(1 + datum.z))
    list(idx = datum.idx, x = datum.x, y = datum.y, seed = datum.seeding, 
        seed.idx = seed.indexing)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

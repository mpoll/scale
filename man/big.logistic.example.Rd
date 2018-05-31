\name{big.logistic.example}
\alias{big.logistic.example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Large Logistic Regression Example
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
big.logistic.example()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    curr.seed <- .Random.seed
    set.seed(1)
    dsz <<- 10^10
    dimen <<- 5
    beta.truth <<- c(1, 1, -1, 2, -2)
    library("msm")
    design.thresh <<- 1
    datum <<- function(datum.idx) {
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
    datum <<- function(datum.idx) {
        old.seed <- .Random.seed
        on.exit({
            .Random.seed <<- old.seed
        })
        seed.indexing <- datum.idx\%/\%10^7
        datum.seeding <- datum.idx\%\%10^7
        set.seed(datum.seeding)
        runif((2 * dimen - 1) * seed.indexing)
        datum.x <- c(1, rtnorm(dimen - 1, mean = 0, sd = 1, lower = -design.thresh, 
            upper = design.thresh))
        datum.z <- exp(sum(beta.truth * datum.x))
        datum.y <- rbinom(1, 1, datum.z/(1 + datum.z))
        list(idx = datum.idx, x = datum.x, y = datum.y, seed = datum.seeding, 
            seed.idx = seed.indexing)
    }
    cycle.size <- 10^7
    cycle.size <- 2147483647
    datum.idx\%/\%cycle.size
    datum.idx\%\%cycle.size
    n.sigma <<- diag(rep(1, dimen)) * dsz^(-1/2)
    n.sigma[1, 1] <<- n.sigma[1, 1]/2
    n.sigma.inv <<- solve(n.sigma)
    beta.star <<- matrix(beta.truth, 1, dimen)
    alpha.cent <<- matrix(c(-0.5848635, -0.3187415, -0.2912184, 
        0.606722, 0.4416645), 1, 5)
    alpha.cent.sq <<- 1.091649
    alpha.p.cent <<- -0.2739236
    phi.cent <<- 0.4088625
    data.extrema()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

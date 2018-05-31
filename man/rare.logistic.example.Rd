\name{rare.logistic.example}
\alias{rare.logistic.example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A variably informative logistic regression example
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rare.logistic.example()
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
    dsz <<- 10^7
    dimen <<- 4
    beta.truth <<- c(0, 2, -2, 2)
    library("msm")
    design.threshold.1 <<- 0.001
    design.threshold.2 <<- 1
    design.thresh <<- c(1, design.threshold.1, rep(design.threshold.2, 
        dimen - 2))
    examp.design <<- cbind(rep(1, dsz), rtnorm(dsz, mean = 0, 
        sd = 1, lower = -design.threshold.1, upper = design.threshold.1), 
        matrix(rtnorm((dimen - 2) * dsz, mean = 0, sd = 1, lower = -design.threshold.2, 
            upper = design.threshold.2), dsz, dimen - 2))
    generate.datum <- function(x) {
        datum.z <- exp(sum(beta.truth * x))
        rbinom(1, 1, datum.z/(1 + datum.z))
    }
    examp.data <<- apply(examp.design, 1, generate.datum)
    sum(examp.data)
    sum(examp.design[, 2])
    initialisation <- logistic.centering()
    beta.star <<- matrix(as.vector(initialisation$center), 1, 
        dimen)
    n.sigma <<- diag(sqrt(as.vector(diag(initialisation$precon))))
    n.sigma.inv <<- solve(n.sigma)
    datum <<- function(datum.idx) {
        datum.x <- examp.design[datum.idx, ]
        datum.y <- examp.data[datum.idx]
        list(idx = datum.idx, x = datum.x, y = datum.y)
    }
    data.init()
    data.extrema()
    .Random.seed <<- curr.seed
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{rare.logistic.example.2}
\alias{rare.logistic.example.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rare.logistic.example.2()
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
    beta.truth <<- c(-20, 2, -2, 2)
    library("msm")
    design.thresh <<- 3
    examp.design <<- cbind(1, matrix(rtnorm((dimen - 1) * dsz, 
        mean = 0, sd = 1, lower = -design.thresh, upper = design.thresh), 
        dsz, dimen - 1))
    generate.datum <- function(x) {
        datum.z <- exp(sum(beta.truth * x))
        rbinom(1, 1, datum.z/(1 + datum.z))
    }
    examp.data <<- apply(examp.design, 1, generate.datum)
    sum(examp.data)
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

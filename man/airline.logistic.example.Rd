\name{airline.logistic.example}
\alias{airline.logistic.example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Airline Data Logistic Regression Example%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
airline.logistic.example()
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
    load("airline_raw.RData")
    dsz <<- 120748239
    dimen <<- 4
    airline <- airline[sample.int(dsz, dsz, replace = FALSE), 
        ]
    examp.design <- matrix(0, dsz, dimen)
    examp.design[, 1] <- 1
    for (j in 2:dimen) {
        examp.design[, j] <- as.numeric(airline[, j])
    }
    examp.design <<- examp.design
    examp.data <<- airline[, 1]
    subsets <- 13
    subset.size <- dsz/subsets
    data.partitions <- 1 + c(0:(subsets - 1)) * subset.size
    subsets.fit <- list()
    for (i in 1:subsets) {
        subset.centering <- logistic.centering(data.partitions[i], 
            subset.size)
        subsets.fit[[i]] <- list(subset = i, idx.start = data.partitions[i], 
            idx.end = data.partitions[i] + subset.size - 1, center = subset.centering$center, 
            precon = subset.centering$precon, precon.inv = subset.centering$precon.inv)
        print(i)
    }
    subset.centers <- matrix(0, subsets, dimen)
    for (i in 1:subsets) {
        subset.centers[i, ] <- subsets.fit[[i]]$center
    }
    beta.star <<- matrix(as.vector(colMeans(subset.centers)), 
        1, dimen)
    precon.inv <- matrix(0, dimen, dimen)
    for (i in 1:subsets) {
        precon.inv <- precon.inv + subsets.fit[[i]]$precon.inv
    }
    precon <- solve(precon.inv)
    n.sigma <<- diag(sqrt(as.vector(diag(precon))))
    n.sigma.inv <<- solve(n.sigma)
    datum <<- function(datum.idx) {
        datum.x <- as.numeric(examp.design[datum.idx, ])
        datum.y <- as.numeric(examp.data[datum.idx])
        list(idx = datum.idx, x = datum.x, y = datum.y)
    }
    subsets <- 13
    subset.size <- dsz/subsets
    data.partitions <- 1 + c(0:(subsets - 1)) * subset.size
    subsets.control <- list()
    for (i in 1:subsets) {
        subset.init <- data.init.2(beta.star, data.partitions[i], 
            subset.size)
        subsets.control[[i]] <- list(subset = i, idx.start = data.partitions[i], 
            idx.end = data.partitions[i] + subset.size - 1, alpha.cent = subset.init$alpha.cent, 
            alpha.p.cent = subset.init$alpha.p.cent)
        print(i)
    }
    subset.alpha.cent <- matrix(0, subsets, dimen)
    subset.alpha.p.cent <- numeric(subsets)
    for (i in 1:subsets) {
        subset.alpha.cent[i, ] <- subsets.control[[i]]$alpha.cent
        subset.alpha.p.cent[i] <- subsets.control[[i]]$alpha.p.cent
    }
    alpha.cent <<- matrix(as.numeric(apply(subset.alpha.cent, 
        2, sum)), nrow = 1)
    alpha.cent.sq <<- (alpha.cent) \%*\% t(alpha.cent)
    alpha.p.cent <<- sum(subset.alpha.p.cent)
    phi.cent <<- (alpha.cent.sq + alpha.p.cent)/2
    subsets <- 13
    subset.size <- dsz/subsets
    data.partitions <- 1 + c(0:(subsets - 1)) * subset.size
    subsets.extreme <- list()
    for (i in 1:subsets) {
        subset.ext <- data.extrema(data.partitions[i], subset.size)
        subsets.extreme[[i]] <- list(subset = i, idx.start = data.partitions[i], 
            idx.end = data.partitions[i] + subset.size - 1, design.max = subset.ext$design.max, 
            design.min = subset.ext$design.min)
        print(i)
    }
    subset.design.min <- subset.design.max <- matrix(0, subsets, 
        dimen)
    for (i in 1:subsets) {
        subset.design.min[i, ] <- subsets.extreme[[i]]$design.min
        subset.design.max[i, ] <- subsets.extreme[[i]]$design.max
    }
    design.min <- as.numeric(apply(subset.design.min, 2, min))
    design.max <- as.numeric(apply(subset.design.max, 2, max))
    data.extrema.2(design.min, design.max)
    .Random.seed <<- curr.seed
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{logistic.centering}
\alias{logistic.centering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute centering value for logistic regression
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logistic.centering(datum.idx.start = 1, length.idxs = dsz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datum.idx.start}{
%%     ~~Describe \code{datum.idx.start} here~~
}
  \item{length.idxs}{
%%     ~~Describe \code{length.idxs} here~~
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
function (datum.idx.start = 1, length.idxs = dsz) 
{
    call.data <- data.subset(datum.idx.start, length.idxs)
    subset.design.mat <- as.matrix(call.data$subset.design[, 
        -1, drop = FALSE])
    subset.data <- call.data$subset.data
    subset.glm <- glm(formula = subset.data ~ subset.design.mat, 
        family = binomial(link = "logit"))
    center <- subset.glm$coef
    precon <- vcov(subset.glm)
    list(center = center, precon = precon, precon.inv = solve(precon), 
        logistic.glm = subset.glm, det.metric = det(precon))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

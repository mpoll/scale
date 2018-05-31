\name{data.subset}
\alias{data.subset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Subset
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.subset(datum.idx.start = 1, length.idxs = dsz)
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
    if (exists("examp.design") == FALSE) {
        call.data <- datum(datum.idx.start)
        subset.design <- t(as.matrix(call.data$x))
        subset.data <- call.data$y
        if (length.idxs > 1) {
            for (i in (datum.idx.start + 1):(datum.idx.start + 
                length.idxs - 1)) {
                call.data <- datum(i)
                subset.design <- rbind(subset.design, call.data$x)
                subset.data <- c(subset.data, call.data$y)
            }
        }
    }
    if (exists("examp.design") == TRUE) {
        subset.design <- examp.design[datum.idx.start:(datum.idx.start + 
            length.idxs - 1), ]
        subset.data <- examp.data[datum.idx.start:(datum.idx.start + 
            length.idxs - 1)]
    }
    list(subset.design = subset.design, subset.data = subset.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

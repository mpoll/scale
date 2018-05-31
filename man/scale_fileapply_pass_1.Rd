\name{scale_fileapply_pass_1}
\alias{scale_fileapply_pass_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_fileapply_pass_1(fnm = "combined_data_1.RData", pathformstart = "scale_precomp_1_", directory = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnm}{
%%     ~~Describe \code{fnm} here~~
}
  \item{pathformstart}{
%%     ~~Describe \code{pathformstart} here~~
}
  \item{directory}{
%%     ~~Describe \code{directory} here~~
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
function (fnm = "combined_data_1.RData", pathformstart = "scale_precomp_1_", 
    directory = getwd()) 
{
    file.names <- dir(path = directory, pattern = pathformstart, 
        full.names = TRUE)
    combined.data_1 <- list()
    for (i in 1:length(file.names)) {
        flist <- NULL
        flist <- ls()
        load(file.names[i], envir = .GlobalEnv)
        run <- as.integer(substring(file.names[i], nchar(directory) + 
            nchar(pathformstart) + 2, nchar(file.names[i]) - 
            nchar(".RData")))
        combined.data_1[[i]] <- list(run = run, subset = subset, 
            idx.start = idx.start, idx.end = idx.end, center = center, 
            precon = precon, precon.inv = precon.inv)
        rm(list = ls()[ls() \%in\% flist == FALSE])
    }
    save(file = fnm, combined.data_1 = combined.data_1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

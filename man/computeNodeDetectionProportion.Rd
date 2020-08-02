\name{computeNodeDetectionProportion}
\alias{computeNodeDetectionProportion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{computeNodeDetectionProportion
%%  ~~function to do ... ~~
}
\description{Takes an input a database table of paths and experimental data of proteomics. The output is the proportion of nodes in a path that were detected by proteomics. It's a useful filter to eliminate paths where the nodes were not detected experimentally.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeNodeDetectionProportion(pdb, sty)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdb}{
%%     ~~Describe \code{pdb} here~~
}
  \item{sty}{
%%     ~~Describe \code{sty} here~~
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
function (pdb, sty)
{
    styDetected <- unique(sty[!is.na(sty$norm), ]$pSite)
    pdb$i2 <- 0
    pdb$i3 <- 0
    pdb$i4 <- 0
    pdb$i5 <- 0
    pdb[pdb$Node2 \%in\% styDetected, ]$i2 <- 1
    pdb[pdb$Node3 \%in\% styDetected & pdb$pathLength \%in\% c(4,
        5, 6), ]$i3 <- 1
    pdb[pdb$Node4 \%in\% styDetected & pdb$pathLength \%in\% c(5,
        6), ]$i4 <- 1
    pdb[pdb$Node5 \%in\% styDetected & pdb$pathLength \%in\% c(6),
        ]$i5 <- 1
    pdb$iSum <- 0
    pdb$iSum <- rowSums(pdb[, c("i2", "i3", "i4", "i5")], na.rm = TRUE)
    pdb$detectedPercent <- pdb$iSum/(pdb$pathLength - 2)
    return(pdb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

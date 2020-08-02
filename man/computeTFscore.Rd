\name{computeTFscore}
\alias{computeTFscore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{computeTFscore
%%  ~~function to do ... ~~
}
\description{The function takes as input two dataframes:
(1) a table database with 2 columns (TranscriptionFactor, Target) both as gene symbols.
(2) an experimental transcriptomics with 2 columns (Gene, norm). The norm is a numeric intensity typically log2 ratios of case / control.
It outputs a Transcription Factor score that is simply the mean expression of all targets.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeTFscore(tdb, rna)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tdb}{
%%     ~~Describe \code{tdb} here~~
}
  \item{rna}{
%%     ~~Describe \code{rna} here~~
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
function (tdb, rna)
{
    tdb2 <- merge(tdb, rna, by.y = "gene", by.x = "Target", all.x = TRUE,
        all.y = FALSE)
    tdb2[tdb2$ControlType \%in\% "Repression", ]$norm <- -1 * tdb2[tdb2$ControlType \%in\%
        "Repression", ]$norm
    t1 <- ddply(tdb2, .(TF), summarize, TFscore = mean(norm,
        na.rm = TRUE))
    t2 <- as.data.table(t1[!is.na(t1$TFscore) & !t1$TFscore ==
        0, ])
    setkey(t2, TF)
    return(t2)
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

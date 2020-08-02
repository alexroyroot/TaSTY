\name{scorePaths}
\alias{scorePaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to find activated paths from receptor kinases to transcription factors ~~
%%  ~~function to find activated paths from receptor kinases to transcription factors ~~
}
\description{
  ~~ The function takes several inputs: (1) a table of paths from receptors to transcription factors (pdb);
 (2) a table of transcription factors and their target genes (tdb)
 (3) experimental data table of phosphoproteomics ratios between 2 experimental conditions (sty) ~~
 (4) experimental data table of proteomics ratios between 2 experimental conditions (pro)
 (5) a constant for the weight of the proteomics data
 (6) a constant for the weight of the phosphoproteomics data
 (7) a constant for the weight of the transcription factor targets
 (8) a cutoff for the percentage of a pathway that must be detected
 (9) the number of permutations to perform
}
\usage{
scorePaths(pdb, tdb, pro,sty, c_pro, c_sty, c_tdb, detectionCutoff, nPerms)
}
%- maybe also 'usage' for other objects documented here.
\details{The function may run very slowly for large datasets with many permutations. Note: the table of paths and the table with transcription factor targets must be wrangled first. The experimental data (proteomics and phosphoproteomics are 3 column dataframes with the column names: gene, norm, norm.sd. The gene column must be valid HGNC symbols. The norm represent the intensity and are typically log2 ratios of two conditions, averaged across biological and techincal replicates. The norm.sd is the standard deviation of these values.

The paths table database (pdb) is a dataframe consisting of 14 columns (Node1, Node2, Node3, Node4, Node5, Node6, pathLength, Gene2, Gene3, Gene4, Gene5, Gene6, TF, pTF))
}
\value{The function returns a data frame with paths, their score, and statistical significance.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Alex Root
%%  ~~who you are~~
}
\note{ This is the main function in the TaSTY package and depends on the data.table package.
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
    print("This is the scorePaths function")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

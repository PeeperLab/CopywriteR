\name{plotCNA}
\alias{plotCNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotCNA segments and plots the copy number files generated with ENCODER.
}
\description{
ENCODER output is analysed by teh segmentation algorithm CBS and plots each sample for each sample the chromosomes seperately. 
}
\usage{
plotCNA(outputFolderENCODER)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outputFolderENCODER}{Path to the output folder of ENCODER}
}
\details{

}
\value{
\item{Plots}{PNG and or PDF figures of the DNA copy number files}
\item{calls.RData}{S4 object with results generated with CGHcall}
}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Guotai Xu, Julian de Ruiter, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F. Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman Submitted for publication. 
}
\author{
%%  ~~who you are~~
Thomas Kuilman (t.kuilman@nki.nl)
Oscar Krijgsman (o.krijgsman@nki.nl)
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
%function (x) 
%{
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

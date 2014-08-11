\name{inputStructure}
\alias{inputStructure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
inputStructure: generates the input files for ENCODER
}
\description{

}
\usage{
ENCODER <- function(inputFile, outputFolder, reference, binSize, bedfile, ncpu){
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{file generated with inputStructure}
  \item{inputFile}{text file with }
  \item{outputFolder}{location of data output (OPTIONAL)}
  \item{reference}{Genome reference reads are mapped to, e.g. hg19 (human), mm10 (mouse).}
  \item{binSize}{Size of the bins for DNA copy number analysis}
  \item{bedfile}{Array design in .bed format. OPTIONAL}
  \item{ncpu}{number of CPUs available for calculation}
}
\details{
For a bin size of 20kb pre-calculated mappability, GC-content and target regions are available within the package. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{inputFile}{input files for ENCODER}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Thomas Kuilman and Oscar Krijgsman
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
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

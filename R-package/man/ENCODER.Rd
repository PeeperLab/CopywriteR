\name{ENCODER}
\alias{ENCODER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ENCODER: off-target sequence reads for DNA copy number detection with whole-exome sequencing
}
\description{

}
\usage{
ENCODER <- function(inputFile){
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{file generated with inputStructure}
}
\details{
DNA copy number profiles directly from whole exome sequencing (WES) based on the captured genomic regions is complicated by the uneven distribution of exons throughout the genome and non-uniform sequence capture. We developed ENCODER (ENhanced COpy number Detection from Exome Reads) which eludes these problems by exploiting the ‘off-target’ sequence reads.
ENCODER counts the reads outside the captured regions and generates a DNA copy number profile with an even distribution and uniform sequence depth over the whole genome. The number of datapoints depends on the number of total reads in the WES bam file. Generally, 

Please note:
Before running ENCODER the function inputStructure should be used to generate the input files for ENCODER. 
For a bin size of 20kb pre-calculated mappability, GC-content and target regions are available within the package. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

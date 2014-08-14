\name{ENCODER}
\alias{ENCODER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ENCODER: off-target sequence reads for DNA copy number detection with whole-exome sequencing
}
\description{

}
\usage{
ENCODER <- function(bamfolder, destinationfolder, referenceFolder, whichControl, captureRegionsBedFile, ncpu){
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfolder}{path to folder with the BAM files}
  \item{destinationfolder}{path to folder where output should be written to}
  \item{referenceFolder}{path to folder with the annotation files generated with preENCODER}
  \item{whichControl}{file that indicates which samples should be used as reference}
  \item{captureRegionsBedFile}{path to the capture design, should have the bed extension}
  \item{ncpu}{number of CPUs available for calculations}
}
\details{
DNA copy number profiles directly from whole exome sequencing (WES) based on the captured genomic regions is complicated by the uneven distribution of exons throughout the genome and non-uniform sequence capture. We developed ENCODER (ENhanced COpy number Detection from Exome Reads) which eludes these problems by exploiting the ‘off-target’ sequence reads.
ENCODER counts the reads outside the captured regions and generates a DNA copy number profile with an even distribution and uniform sequence depth over the whole genome. The number of datapoints depends on the number of total reads in the WES bam file. Generally, a 20kb bin size works for most WES samples. 
Before running ENCODER the function inputStructure should be used to generate the input files for ENCODER. 




}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{log.txt}{log file of ENCODER}
\item{BamBaiMacsFiles}{Folder with the files made by ENCODER needed to run the algorithm}
\item{CNAprofiles}{Folder with the output generate by ENCODER}
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
#ENCODER <- function(bamfolder, destinationfolder, referenceFolder, whichControl, captureRegionsBedFile, ncpu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{preENCODER}
\alias{preENCODER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
preENCODER is used to generate annotation files for the desired bin size
}
\description{
In ENCODER genome-wide annotation files with mappability, GC_content and blacklists for 1000kb bins are providedmin binary .rda files. The function preENCODER reads the provided data files and generates annotation files for the desired bin size. The output folder will contain the files needed to run the main function ENCODER. 
For each reference genome and binSize this function has to be executed once. 
}
\usage{
preENCODER(MAPA_GC_location, outputFolder, binSize, reference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAPA_GC_location}{path to the 1kb annotation files}
  \item{outputFolder}{path to the output folder}
  \item{binSize}{desired bin length for which annotation files should be generated}
  \item{reference}{reference, e.g. hg19 or mm10}
}
\details{
Currently annotation files can be generated for both human (hg19) and mouse (mm10) genomes. Each 
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
Oscar krijgsman (o.krijgsman@nki.nl)
Thomas Kuilman (t.kuilman@nki.nl)
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

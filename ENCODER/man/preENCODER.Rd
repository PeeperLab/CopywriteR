\name{preENCODER}
\alias{preENCODER}
\title{
preENCODER is used to generate annotation files for the desired bin size
}
\description{
In ENCODER genome-wide annotation files with mappability, GC_content and blacklists for 1000kb bins are provided as binary .rda files. The function preENCODER reads the provided data files and generates annotation files for the desired bin size. The output folder will contain the files needed to run the main function ENCODER. 
For each reference genome and binSize this function has to be executed once. 
}
\usage{
preENCODER(MAPA_GC_location, outputFolder, binSize, reference)
}
\arguments{
  \item{MAPA_GC_location}{path to the 1kb annotation files}
  \item{outputFolder}{path to the output folder}
  \item{binSize}{desired bin length for which annotation files should be generated}
  \item{reference}{reference, e.g. hg19 or mm10}
}
\details{
Currently annotation files can only be generated for  human (hg19). 
We are currently implementing mm10 and will be available in the next release. 

Annotation files can only be generated with binSize that are a multiplication of 1000bp (1kb). 
}
\value{
\item{folder}{folder with generated annotation files with desired binSize}
\item{input.Rdata}{R object with all generated data}

}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Guotai Xu, Julian de Ruiter, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F. Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman Submitted for publication. 
}
\author{
Oscar krijgsman (o.krijgsman@nki.nl)
Thomas Kuilman (t.kuilman@nki.nl)
}

\examples{
## 
\dontrun{preENCODER("./annotation/", "./", 20000, "hg19")}
}
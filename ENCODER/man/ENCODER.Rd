\name{ENCODER}
\alias{ENCODER}
\title{
ENCODER: off-target sequence reads for DNA copy number detection with whole-exome sequencing
}
\description{
Generates DNA copy number profiles from whole exome sequencing bam files
}
\usage{
ENCODER(bamFolder, destinationFolder, referenceFolder, whichControl, ncpu, captureRegionsBedFile)
}

\arguments{
  \item{bamFolder}{path to folder with the BAM files}
  \item{destinationFolder}{path to folder where output should be written to}
  \item{referenceFolder}{path to folder with the annotation files generated with \code{\link{preENCODER}} }
  \item{whichControl}{file that indicates which samples should be used as reference}
  \item{captureRegionsBedFile}{path to the capture design, should have the bed extension}
  \item{ncpu}{number of CPUs available for calculations}
}
\details{
DNA copy number profiles directly from whole exome sequencing (WES) based on the captured genomic regions is complicated by the uneven distribution of exons throughout the genome and non-uniform sequence capture. We developed ENCODER (ENhanced COpy number Detection from Exome Reads) which eludes these problems by exploiting the off-target sequence reads.
ENCODER counts the reads outside the captured regions and generates a DNA copy number profile with an even distribution and uniform sequence depth over the whole genome. The number of datapoints depends on the number of total reads in the WES bam file. Generally, a 20kb bin size works for most WES samples. 
Before running ENCODER the function inputStructure should be used to generate the input files for ENCODER. 

}
\value{
  \item{log.txt}{log file of ENCODER}
  \item{BamBaiMacsFiles}{Folder with the files made by ENCODER needed to run the algorithm}
  \item{CNAprofiles}{Folder with the output generate by ENCODER}
}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Guotai Xu, Julian de Ruiter, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F. Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman Submitted for publication. 
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
Oscar Krijgsman (o.krijgsman@nki.nl)
}

\examples{
## Generate DNA copy number data from whole exome sequence data with annotation as generated with \code{\link{preENCODER}}.
  \dontrun{ENCODER("/", "/output", "/hg19_20kb", c(2,2,2), capturefile.bed, 4)}
}

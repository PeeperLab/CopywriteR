\name{plotCNA}
\alias{plotCNA}
\title{
ENCODER: ENhanced COpy number Detection from Exome Reads
}
\description{
plotCNA analyses ENCODER output using the segmentation algorithm CBS and creates
both whole-genome and per-chromosome plots.
}
\usage{
plotCNA(destination.folder, smoothed = TRUE, sample.plot, y.min, y.max, ...)
}
\arguments{
\item{destination.folder}{
the path to the output folder of ENCODER.
}
\item{smoothed}{
logical that determines whether DNAcopy smoothens copy number information before
segmentation. Defaults to TRUE.
}
\item{sample.plot}{
optional; a data.frame or matrix that contains the locations of sample and
control BAM files. The column names should be "samples" and "controls",
respectively. The file locations in a row represent a sample and its
corresponding control, which is used for plotting of relative copy number
values. Every sample is plotted both without a reference and with a reference.
Defaults to sample.controls as used for running ENCODER.
}
\item{y.min}{
determines lower boundary of the plotting range.
}
\item{y.max}{
determines upper boundary of the plotting range.
}
\item{...}{
takes additional arguments that will be passed on to DNAcopy.
}
}
\value{
\item{segment.Rdata}{
an R object of the DNAcopy class, which contains the segmentation values for
all analysed samples.
}
\item{plots}{
a folder that contains all copy number profiles in portable document file
(pdf) file format. Both profiles per chromosome, as well as a genome-wide
profile are provided for every sample.
}
}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas
Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli,
Marlous Hoogstraat, Ekaterina Nevedomskaya, Guotai Xu, Julian de Ruiter,
Martijn P. Lolkema, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F.
Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman. Submitted for
publication. 
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
}
\seealso{
\code{\link{ENCODER}}
}
\examples{
\dontrun{
setwd("/PATH/TO/BAMFILES/")
samples <- list.files(pattern = ".bam$", full.names = TRUE)
controls <- samples[rep(1, length(samples))] ## Use the first .bam file as a reference for plotting every sample
sample.plot <- data.frame(samples, controls)

plotCNA("./PATH/TO/DESTINATIONFOLDER/", sample.plot = sample.plot, y.min = -3,
        y.max = 3)
}
}

  
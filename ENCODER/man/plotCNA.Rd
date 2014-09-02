\name{plotCNA}
\alias{plotCNA}
\title{
plotCNA segments and plots the copy number files generated with ENCODER.
}
\description{
ENCODER output is analysed by the segmentation algorithm CBS and creates both whole-genome and per-chromosome plots.
}
\usage{
plotCNA(destinationFolder, set.nchrom)
}
\arguments{
  \item{destinationFolder}{Path to the output folder of ENCODER}
  \item{set.nchrom}{Number of chromosomes that are included in analysis. Use this option if CGHcall/plotCNA fail due to insufficient data points on the Y-chromosome.}
}

\value{
\item{Plots}{PNG and or PDF figures of the DNA copy number files}
\item{calls.RData}{S4 object with results generated with CGHcall}
}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Guotai Xu, Julian de Ruiter, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F. Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman. Submitted for publication. 
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
Oscar Krijgsman (o.krijgsman@nki.nl)
}

\seealso{
\code{\link{ENCODER}}
}
\examples{
## plot DNA copy number profiles from \code{\link{ENCODER}} output.
  \dontrun{plotCNA("./PATH_TO_DESTINATIONFOLDER/")}
  \dontrun{plotCNA("./PATH_TO_DESTINATIONFOLDER/", set.nchrom = 23)} ## Use the total chromosome count minus 1 (i.e., 23 for the human genome) when CGHcall fails
  }

  
\name{plotCNA}
\alias{plotCNA}
\title{
plotCNA segments and plots the copy number files generated with ENCODER.
}
\description{
ENCODER output is analysed by the segmentation algorithm CBS and plots each sample for each sample the chromosomes seperately. 
}
\usage{
plotCNA(destinationfolder)
}
\arguments{
  \item{destinationfolder}{Path to the output folder of ENCODER}
}

\value{
\item{Plots}{PNG and or PDF figures of the DNA copy number files}
\item{calls.RData}{S4 object with results generated with CGHcall}
}
\references{
ENCODER: DNA copy number detection from off-target sequence data. Thomas Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Guotai Xu, Julian de Ruiter, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F. Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman Submitted for publication. 
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
  \dontrun{plotCNA(./destinationFolder/)}
  }

  
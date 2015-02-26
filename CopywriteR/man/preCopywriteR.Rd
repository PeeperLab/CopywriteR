\name{preCopywriteR}
\alias{preCopywriteR}
\title{
CopywriteR: DNA copy number detection from off-target sequence data
}
\description{
preCopywriteR is used to generate helper files (blacklist, bin region, GC-content,
and mappability .bed files) for the desired bin size from pre-assembled 1 kb 
helper files. These binary 1 kb helper files are available at
\url{https://github.com/PeeperLab/CopywriteR} for hg19, mm9 and mm10 reference
genomes.
}
\usage{
preCopywriteR(output.folder, bin.size, ref.genome)
}
\arguments{
\item{output.folder}{
the path to the output folder.
}
\item{bin.size}{
desired bin length (in bp) for which helper files should be generated.
}
\item{ref.genome}{
the reference genome for which helper files should be generated, e.g. hg19, mm9
or mm10.
}
}
\details{
Currently helper files can only be generated for  human (hg19). We are
currently implementing mm10 and will be available in the next release.

Helper files can only be generated for a bin.size that is a multiple of 1000 bp. 
}
\value{
\item{blacklist.rda}{
an R data file containing a GRanges object with blacklisted regions of known
CNVs.
}
\item{GC_mappability.rda}{
an R data file containing a GRanges object with the mappability and GC-content
for bins of the specified size.
}
}
\references{
CopywriteR: DNA copy number detection from off-target sequence data. Thomas
Kuilman, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli,
Marlous Hoogstraat, Ekaterina Nevedomskaya, Guotai Xu, Julian de Ruiter,
Martijn P. Lolkema, Bauke Ylstra, Jos Jonkers, Sven Rottenberg, Lodewyk F.
Wessels, David J. Adams, Daniel S. Peeper, Oscar Krijgsman. Submitted for
publication. 
}
\author{
Oscar krijgsman (o.krijgsman@nki.nl)\cr
Thomas Kuilman (t.kuilman@nki.nl)
}
\examples{
\dontrun{
preCopywriteR("/PATH/TO/HG19_1KB_FOLDER/", "./", 20000, "hg19")
}
}
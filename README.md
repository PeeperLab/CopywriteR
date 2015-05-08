# CopywriteR
#### (formerly known as [ENCODER](https://github.com/PeeperLab/ENCODER))

Current methods for detection of copy number variants and aberrations (CNV and
CNA) from targeted sequencing data are based on the depth of coverage of
captured exons. Accurate CNA determination is complicated by uneven genomic
distribution and non-uniform capture efficiency of targeted exons. Here we
present CopywriteR which eludes these problems by exploiting ‘off-target’
sequence reads. CopywriteR allows for extracting uniformly distributed copy
number information, can be used without reference and can be applied to
sequencing data obtained from various techniques including chromatin
immunoprecipitation and target enrichment on small gene panels. CopywriteR
outperforms existing methods and constitutes a widely applicable alternative to
available tools.

CopywriteR has been described in
[Kuilman et al., 2015](http://genomebiology.com/2015/16/1/49/abstract). The
analysis described in this publication were performed using the older version
[V1.3](https://github.com/PeeperLab/CopywriteR/releases/tag/V1.3).

## Bioconductor

We are happy to announce that the CopywriteR package has been accepted in
Bioconductor and has been released as of 17 April 2015. All the Bioconductor
packages in the newest release have a dependency for the newest version of R
(version 3.2) and for Bioconductor version 3.1. If you have these installed,
CopywriteR can be installed as follows:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("CopywriteR")

We have developed CopywriteR with R 3.1 and Bioconductor 3.0, and therefore know
that it works fine with these version too. In case you are running these,
you can still install CopywriteR as is explained below.

## Installation (not via Bioconductor)

CopywriteR can be installed without Bioconductor, which is useful when you have
R 3.1 and Bioconductor 3.0 installed. Installation should be performed as
follows:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("matrixStats", "gtools", "data.table", "S4Vectors", "chipseq",
                 "IRanges", "Rsamtools", "DNAcopy", "GenomicAlignments",
                 "GenomicRanges", "GenomeInfoDb", "BiocParallel",
                 "futile.logger"))

In addition, CopywriteR requires the CopyhelpeR package, which can be downloaded
as tarball (.tar.gz-file) from the
[CopyhelpeR releases webpage](https://github.com/PeeperLab/CopyhelpeR/releases).
Subsequently, the CopyhelpeR package can be installed from the command line
using the following command:

    $ R CMD INSTALL CopyhelpeR*.tar.gz

As the last step in the installation process, the latest CopywriteR package can
be downloaded from the
[CopywriteR releases webpage](https://github.com/PeeperLab/CopywriteR/releases)
and installed using the following command:

    $ R CMD INSTALL CopywriteR*.tar.gz

Now you are all set to start your analysis.

## CopywriteR usage:

Load the CopywriteR package in R using:

    > library("CopywriteR")

CopywriteR contains three main functions:

`preCopywriteR` will generate a GRanges object containing mappability and
GC-content information, and one containing 'blacklisted' regions that contain
are subject to copy number variation. These 'helper' files can be created for
any specified bin size that is a multiple of 1000 bp, and for any of the
available reference genomes (hg18, hg19, hg38, mm9 and mm10). The helper files
can be re-used and need to be created only once for every combination of
reference genome and bin size. `preCopywriteR` uses information stored in
pre-assembled 1kb bin mappability and GC-content GRanges objects to create the
custom bin size helper files. These objects are stored in the CopyhelpeR
annotation package.

preCopywriteR can be run as follows:

    > preCopywriteR(output.folder, bin.size, ref.genome, prefix = "")

`CopywriteR` will generate separate tables with compensated read counts and
log2-transformed normalized compensated read counts after correction for
GC-content, mappability and upon removal of blacklisted regions.

    > CopywriteR(sample.control, destination.folder, reference.folder,
                 bp.param, capture.regions.file,
                 keep.intermediary.files = FALSE)

`plotCNA` performs segmentation using the DNAcopy Bioconductor package, and
plotting of copy number profiles.

    > plotCNA(destinationFolder, set.nchrom)

For more details please refer to the CopywriteR vignette. Alternatively, one of
the following commands can be used to show help files for the corresponding
function:

    > ?preCopywriteR
    > ?CopywriteR
    > ?plotCNA

## Quick start

A typical analysis using CopywriteR could be as follows. First, CopywriteR needs
to be loaded:

    > library(CopywriteR)

Then, preCopywriteR can be run using the command:

    > preCopywriteR(output.folder = file.path("./path/to/output/folder"),
                    bin.size = 20000,
                    ref.genome = "mm10",
                    prefix = "")

Next, we need to specify the settings for parallel computing. We have
implemented use of the
[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
package, which supports several different types of environments. For every
environment, a BiocParallelParam can be specified that defines how parallel
computation is executed. Below, we use a SnowParam instance of
BiocParallelParam, which is based on the
[snow](http://cran.r-project.org/web/packages/snow/index.html) package. Please
refer to the
[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
package for more information. A SnowParam using 12 CPUs can be defined as
follows:

    > bp.param <- SnowParam(workers = 12, type = "SOCK")

Next, we need to specify which samples and controls correspond to each other
using the sample.control variable. **For the CopywriteR function, controls are
specified as those samples that will be used to identify which regions are
'peaks' and contain on-target reads.** This information will then be used to
remove on-target reads in the corresponding sample. We specify the
sample.control variable as follows:

    > path <- "./path/to/bam"
    > samples <- list.files(path = path, pattern = "bam$", full.names = TRUE)
    > controls <- samples[c(1:6, 1:6)]
    > sample.control <- data.frame(samples, controls)

This might result in the following variable:

    > sample.control
          samples   controls
    1  ./C003.bam ./C003.bam
    2  ./C016.bam ./C016.bam
    3  ./C024.bam ./C024.bam
    4  ./C037.bam ./C037.bam
    5  ./C049.bam ./C049.bam
    6  ./C055.bam ./C055.bam
    7  ./M003.bam ./C003.bam
    8  ./M016.bam ./C016.bam
    9  ./M024.bam ./C024.bam
    10 ./M037.bam ./C037.bam
    11 ./M049.bam ./C049.bam
    12 ./M055.bam ./C055.bam

Sequence data starting with 'M' could for instance be from a tumor sample, and
the corresponding 'C' data set would from a matched germline sample. **Please
note that any sample that is to be used by the downstream plotCNA function needs
to be analyzed by the CopywriteR function.** Therefore, by including:

    1  ./C003.bam ./C003.bam
    2  ./C016.bam ./C016.bam
    3  ./C024.bam ./C024.bam
    4  ./C037.bam ./C037.bam
    5  ./C049.bam ./C049.bam
    6  ./C055.bam ./C055.bam

we make sure that in the plotCNA function we can analyze the tumor samples
relative to the corresponding germline samples. We recommend identifying
on-target and off-target regions based on a germline sample if possible, as this
would avoid identifying highly amplified genomic regions in tumor cells as
on-target regions. Nevertheless, we have observed that this effect is negligible
in practice, and that CopywriteR analysis without a reference is still highly
accurate. Please refer to
[Kuilman et al., 2015](http://genomebiology.com/2015/16/1/49/abstract) for more
details on how CopywriteR extracts copy number profiles from targeted sequencing.

We are now set for running CopywriteR:

    > CopywriteR(sample.control = sample.control,
                 destination.folder = file.path("./path/to/destination/folder"),
                 reference.folder = file.path("./path/to/reference/folder", "mm10_20kb"),
                 bp.param = bp.param)

Finally, we can segment the data using
[DNAcopy](http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html)
and plot the result using the plotCNA function:

    > plotCNA(destination.folder = file.path("/path/to/destination/folder"))

By default, data will be analyzed and plotted according to the matching of
samples and controls in the sample.control variable (specified above). Every
sample in the samples column of the sample.control variable will be analyzed
without a reference, and with the corresponding control as a reference.
Optionally, the sample.plot argument can be used to control analysis and
plotting by plotCNA. Please refer to the manual for more information.

## Vignette code

The following packages can be optionally installed to allow running the
vignette code:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("BiocStyle", "snow"))

In addition, the SCLCBam experiment package is also required. The latest
SCLCBam package can be downloaded from the
[SCLCBam releases webpage](https://github.com/PeeperLab/SCLCBam/releases)
and installed using the following command:

    $ R CMD INSTALL SCLCBam*.tar.gz

## Troubleshooting
In general, we advise to update to the latest versions of CopywriteR and
CopyhelpeR in case of errors. If the problems persist, please refer below for
troubleshooting purposes.

#### Mm10
We have come to realise that there was an bug in CopyhelpeR version 1.0.0 which
leads to an error when creating mm10 helper files using the preCopywriteR
function. If you need to use CopywriteR for mm10-based sequence data, please
update to version 1.0.1 available at
[GitHub](https://github.com/PeeperLab/CopyhelpeR/releases) or via
[bioconductor](http://bioconductor.org/packages/release/data/experiment/html/CopyhelpeR.html)
(available with the next build).

#### Installation
One of the dependencies of CopywriteR (the chipseq package) requires
Bioconductor 3.0. If installation fails, please check whether you are running
the correct version of Bioconductor and whether all dependencies have been
installed. When using the BiocInstaller

#### Unique naming of .bam files
CopywriteR uses the (base)name of a .bam file as an identifier for that file.
Therefore, all the names of .bam files should be unique, and something like

    > sample.control
                             samples                        controls
    1  path/sample1/aligned_read.bam   path/sample2/aligned_read.bam

would result in an error. This can be solved by renaming the .bam files in such
a way that all names become unique.

#### All .bam files should be processed by CopywriteR
Any sample that is used as a sample or a reference for analysis and plotting in
the plotCNA function needs prior analysis as a sample in the CopywriteR
function. As an example, if one would like to analyze and plot tumor1.bam
relative to matched.normal1.bam in plotCNA, the sample.control variable should
contain both of the two rows below:

    > sample.control
                           samples                    controls
    1  path/to/tumor1.bam          path/to/matched.normal1.bam
    2  path/to/matched.normal1.bam path/to/matched.normal1.bam

If not all samples needed by plotCNA have been analyzed by CopywriteR, the
following error message will be displayed: "One of the samples in [LIST OF .BAM
FILES] refers to a BAM file that has not been processed in CopywriteR. Please
make sure that you have provided the correct input files or re-run CopywriteR
accordingly."

#### Number of off-target reads
CopywriteR relies on the off-target sequence read count, which in our hands is
roughly 10% of the total amount of reads. You can find the relevant numbers in
the log-file. Since the sequence reads are normally not required, a number of
pipelines filter for or deplete sequence reads that are off-target. We recommend
using the unprocessed .bam files to prevent discarding data for CopywriteR.

In addition to the above, we have never tested CopywriteR on sequence data upon
PCR-based enrichment strategies. However, we believe it is likely that there
will be insufficient off-target reads to allow high-quality copy number data
using CopywriteR.

#### Chromosome names

CopywriteR by default assumes that the chromosome names in .bam files are "1",
"2", ... "X", and "Y". These chromosome names are incorporated in the
mappability and GC-content GRanges object, as well as in the blacklist GRanges
object created by preCopywriteR. preCopywriteR has an optional prefix argument
that can be used to created helper files with a chromosome notation using a
prefix (for instance "chr" for "chr1", "chr2", ... notation).

Moreover, CopywriteR assumes that all .bam files that are analyzed together are
mapped to the same reference genome. If this is not the case, please re-run
multiple CopywriteR analyses with every analysis using .bam files mapped to the
same reference genome.

#### Plotting

In order to allow plotting, the .bam files need to be present in the same folder
as when the CopywriteR analysis was performed.

#### Sorting of .bam files

The Rsamtools package requires .bam files to be sorted. If the .bam files that
are to be analyzed are not sorted, CopywriteR will provide an error message.
Unsorted .bam files can then be sorted using samtools as follows:

    $ samtools sort sample.bam sample.sorted.bam

#### Speed

We have noted a ~1.5-fold reduction in total calculation time with the
development version of R (3.2) and Bioconductor (3.1). These will soon be
officially released, and in case the calculation time becomes a critical issue
we recommend upgrading to these versions of R and Bioconductor.

## Contact

We have tried to make the CopywriteR code readable and its use as easy as
possible. If any questions arise regarding the package, or if you want to report
any bugs, please do not hesitate and contact:

- [Thomas Kuilman](mailto:t.kuilman@nki.nl)
- [Oscar Krijgsman](mailto:o.krijgsman@nki.nl)

Thomas and Oscar are working in the laboratory of Prof. Dr. Daniel S. Peeper.

- [Lab website](http://research.nki.nl/peeperlab/)

## Cited by

The CopywriteR tool has been cited and referenced by:

- Hoogstraat, M., Gadellaa-van Hooijdonk, C. G., Ubink, I., Besselink, N. J. M.,
  Pieterse, M., Veldhuis, W., van Stralen, M., Meijer, E. F. J., Willems, S. M.,
  Hadders, M. A., Kuilman, T., Krijgsman, O., Peeper, D. S., Koudijs, M. J.,
  Cuppen, E., Voest, E. E. and Lolkema, M. P. (2015), Detailed imaging and
  genetic analysis reveal a secondary BRAFL505H resistance mutation and
  extensive intrapatient heterogeneity in metastatic BRAF mutant melanoma
  patients treated with vemurafenib. Pigment Cell & Melanoma Research.
  doi: 10.1111/pcmr.12347.
  [[Pubmed]](http://www.ncbi.nlm.nih.gov/pubmed/25515853)
- Schouten, P.C., Grigoriadis, A., Kuilman, T., Mirza, H., Watkins, J.A., Cooke,
  S.A., van Dyk, E., Severson, T.M., Rueda, O.M., Hoogstraat, M., Verhagen, C.,
	Natrajan, R., Chin, S.F., Lips, E.H., Kruizinga, J., Velds, A., Nieuwland, M.,
	Kerkhoven, R.M., Krijgsman, O., Vens, C., Peeper, D., Nederlof, P.M., Caldas,
	C., Tutt, A.N., Wessels, L.F. and Linn, S.C. (2015), Robust BRCA1-like
	classification of copy number profiles of samples repeated across different
	datasets and platforms. Molecular Oncology. doi: 10.1016/j.molonc.2015.03.002.
	[[Pubmed]](http://www.ncbi.nlm.nih.gov/pubmed/25825120)
- [Omicstools](http://omictools.com/copywriter-s8589.html)

## Reported bugs

We have received feedback that CopywriteR throws the following error in a
platform-dependent manner:

Error in log2.read.counts[selection, , drop = FALSE] :
  (subscript) logical subscript too long

We are working on removing this error, but this might take a bit longer as we
cannot reproduce this error on our platform. You can stay updated by 'watching'
the CopywriteR package at the
[watchers](https://github.com/PeeperLab/CopywriteR/watchers) page.

## Changes and additions we are currently working on

- [ ] Include check for empty peaks.bed files when calculating overlap
- [ ] Remove plotting of self vs self
- [ ] Remove blacklisted regions by excluding regions and not excluding bins
- [ ] Change CopywriteR to allow custom bins?
- [x] Compile into bioConductor package (implemented in source code)
- [x] Reset the working directories after analysis
- [x] Fix bug related to bin-spanning peaks (implemented in source code)
- [x] Fix bug related to analysis of non-indexed .bam files (implemented in source code)
- [x] Fix bug related to testing sorting of .bam files (implemented in source code)
- [x] Fix bug affecting read count statistics on single-sample analysis (implemented in source code)
- [x] Remove warning messages small bams _chr (implemented in source code)
- [x] Fixed a minor bug affecting analysis of >100 simultaneous samples (implemented in source code)
- [x] Change plotCNA (double line in plots of RB1705) (implemented in source code)
- [x] Improve error message in plotCNA when CopywriteR analysis is missing for some samples (implemented in source code)
- [x] Clean up code (implemented in source code)
- [x] Use futile.logger for logging of messages (implemented in source code)
- [x] Increase speed of paired / single end testing (implemented in source code)
- [x] Make naming consistent (mappability, chromosome etc) (implemented in source code)
- [x] Solve problem with 0-based ranges (implemented in source code)
- [x] Implement GenomicRanges::makeGRangesFromDataFrame() (implemented in source code)
- [x] Add version of CopywriteR to log (implemented in source code)
- [x] Address all NOTES when building package (implemented in source code)
- [x] Check problem with not using all bams as samples (RB; implemented in source code)
- [x] Include helper file folder in log-file (implemented in source code)
- [x] Provide output naturally sorted (implemented in source code)
- [x] Solve issue with data RB (implemented in source code)
- [x] Remove potential to overwrite existing folders and files (implemented in source code)
- [x] Change handling of NAs in plotCNA function (implemented in source code)
- [x] Make names consistent and apply Rlint to code (implemented in source code)
- [x] Use BiocParallel for parallel computing (implemented in source code)
- [x] Check why sex chromosomes are not plotted (solved in source code)
- [x] Add track line to .igv output file (#track viewLimits=-3:3 graphType=heatmap color=255,0,0) (implemented in source code)
- [x] Check log.txt output file (implemented in source code)
- [x] Check that all chromosomes have same length (implemented in source code)
- [x] Remove possibility for error when allowing more cpus than samples (implemented in source code)
- [x] Check for platform-compatibility when handling sample.files variable (implemented in source code)
- [x] Check potential plotting error when running on custom bin files (implemented in source code)
- [x] Remove dependency for MACS 1.4 and use chipseq package instead (implemented in source code)
- [x] Change -Inf in log2-table to small value for visualization in IGV (implemented in source code)
- [x] Provide the amount of sequence reads used by CopywriteR (implemented in source code)
- [x] Check rm() function if no capture region file is specified (implemented in source code)
- [x] Provide the option to keep intermediary files; remove by default (implemented in source code)
- [x] Check writing permissions on bam.folder (implemented in source code)
- [x] Include check for presence .bai files (implemented in source code)
- [x] Use file.paths for all file paths (implemented in source code)
- [x] Change output format for raw data from .txt to .igv (implemented in source code)
- [x] Remove bug in determining the total chromosome number with missing sex chromosomes (implemented in source code)
- [x] Remove dependency for CGHcall and use only DNAcopy for segmentation (implemented in source code)
- [x] Provide plotting function for DNAcopy output (implemented in source code)
- [x] Make keeping intermediate .bam and MACS files optional (implemented in source code)
- [x] Change input structure to allow more flexibility (implemented in source code)
- [x] Clean up code (implemented in source code)
- [x] Remove dependency for bedtools (replaced by GenomicRanges in source code)
- [x] Remove dependency for Samtools (replaced by Rsamtools in source code)
- [x] Include a check for integrity of BAM files (implemented in source code)
- [x] Change output format to igv2.0 (implemented in source code)
- [x] Remove bug in determining chromosome name prefixes (fixed in source code)
- [x] Debug error when using 1kb bin size
- [x] Support alternative chromosome names (i.e., "chr1" instead of "1")
- [x] Extract binSize from bins.bed file
- [x] Support for relative path names
- [x] Remove requirement for a trailing `/` in folder path names
- [x] Add 1kb mapability and GC-content files for mouse mm9 and mm10 genomes
- [x] Make captureRegionsBedFile optional
- [x] Allow processing of single-end sequences
- [x] Increase speed for generating bins in `preCopywriteR`